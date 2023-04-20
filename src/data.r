library(purrr)
library(stringr)
library(ms.sev)

source("src/edmus.r")
source("src/biobank.r")
source("src/hla.r")


patients <- edmus_personal |>
    select(
        patient_id,
        edmus_local_id = local_identifier,
        nhc = other_identifier,
        gender,
        last_clinical_follow_up,
        date_of_birth,
        wait_and_see,
        starts_with("irreversible_dss_") & !ends_with("_unknown_date")
    ) |>
    rename_with(~ str_replace(.x, "^irreversible_dss_", "date_of_iedss_")) |>
    left_join(
        edmus_diagnosis |>
            select(patient_id, ms_onset, disease_course, progression_onset) |>
            mutate(disease_course = case_match(
                disease_course, "RR" ~ "RR",
                "SP-R" ~ "SP", "SP-NR" ~ "SP",
                "PP-R" ~ "PP", "PP-NR" ~ "PP"
            )),
        by = "patient_id"
    ) |>
    filter(
        !wait_and_see,
        disease_course %in% c("RR", "SP")
    ) |>
    select(-wait_and_see) |>
    mutate(
        disease_duration = last_clinical_follow_up - ms_onset,
    ) |>
    inner_join(
        hla_data |>
            drop_na(codi_edm) |>
            select(codi_edm, starts_with("tipatge_hla_drb1_")) |>
            rename_with(~ str_replace(.x, "^tipatge_hla_drb1_", "hla_drb1_")),
        by = c(edmus_local_id = "codi_edm"), multiple = "first"
    )

rr_relapse_counts <- patients |>
    left_join(
        edmus_episodes |> select(patient_id, date),
        by = "patient_id", relationship = "many-to-many"
    ) |>
    filter(is.na(progression_onset) | date < progression_onset) |>
    group_by(patient_id) |>
    summarize(relapse_count_until_progression = n())

patients <- patients |>
    left_join(rr_relapse_counts, by = "patient_id") |>
    mutate(
        relapse_count_until_progression = replace_na(
            relapse_count_until_progression, 0
        ),
        duration_of_rr_phase = if_else(
            disease_course |> str_starts("SP"),
            progression_onset - ms_onset,
            last_clinical_follow_up - ms_onset
        ),
        arr_during_rr_phase = if_else(
            duration_of_rr_phase >= dyears(1),
            relapse_count_until_progression / (duration_of_rr_phase / dyears(1)),
            NA_real_
        )
    )

clinical_iedss <- edmus_clinical |>
    semi_join(patients, by = "patient_id") |>
    select(patient_id, date, edss_entered) |>
    group_by(patient_id) |>
    filter((max(date) - date) >= dmonths(6)) |>
    arrange(desc(date), .by_group = TRUE) %>%
    mutate(
        iedss = accumulate(edss_entered, \(x, y) {
            if (all(!is.finite(c(x, y)))) {
                return(NA)
            }
            min(x, y, na.rm = TRUE)
        })
    ) |>
    ungroup()

treatment_times <- edmus_trt_dm |>
    mutate(time_on_treatment = case_match(
        status,
        "Stopped" ~ end_date - onset_date,
        "Ongoing" ~ latest_date - onset_date
    )) |>
    group_by(patient_id) |>
    summarize(
        time_on_high_efficacy_treatment = sum(if_else(
            inn %in% edmus_trt_dm_high_efficacy_inn, time_on_treatment, ddays(0)
        )),
        time_on_moderate_efficacy_treatment = sum(if_else(
            inn %in% edmus_trt_dm_moderate_efficacy_inn, time_on_treatment, ddays(0)
        )),
        time_on_other_treatments = sum(if_else(
            !inn %in% c(edmus_trt_dm_high_efficacy_inn, edmus_trt_dm_moderate_efficacy_inn),
            time_on_treatment, ddays(0)
        ))
    )

iedss_global <- clinical_iedss |>
    slice_max(date, by = "patient_id", n = 1, with_ties = FALSE) |>
    select(patient_id, date, iedss_global = iedss)

iedss_rr_phase <- clinical_iedss |>
    left_join(patients, by = "patient_id") |>
    filter(disease_course == "RR" | date <= progression_onset) |>
    slice_max(date, by = "patient_id", n = 1, with_ties = FALSE) |>
    select(patient_id, date, iedss_rr_phase = iedss)

mssev_params_global <- iedss_global |>
    left_join(patients, by = "patient_id") |>
    transmute(
        patient_id,
        edss = iedss_global,
        dd = (date - ms_onset) / dyears(1),
        ageatedss = (date - date_of_birth) / dyears(1)
    )

mssev_params_rr_phase <- iedss_rr_phase |>
    left_join(patients, by = "patient_id") |>
    transmute(
        patient_id,
        edss = iedss_rr_phase,
        dd = (date - ms_onset) / dyears(1),
        ageatedss = (date - date_of_birth) / dyears(1)
    )

armss_global <- ms_sev(mssev_params_global, type = "global_armss")$data |>
    select(patient_id, armss_global = "gARMSS")

armss_rr_phase <- ms_sev(mssev_params_rr_phase, type = "global_armss")$data |>
    select(patient_id, armss_rr_phase = "gARMSS")

msss_global <- ms_sev(mssev_params_global, type = "global_msss", omsss = TRUE)$data |>
    filter(dd >= 1) |>
    select(patient_id, msss_global = "oGMSSS")

msss_rr_phase <- ms_sev(mssev_params_rr_phase, type = "global_msss", omsss = TRUE)$data |>
    filter(dd >= 1) |>
    select(patient_id, msss_rr_phase = "oGMSSS")

severity_scores <- armss_global |>
    left_join(iedss_global |> select(-date), by = "patient_id") |>
    left_join(iedss_rr_phase |> select(-date), by = "patient_id") |>
    left_join(armss_rr_phase, by = "patient_id") |>
    left_join(msss_global, by = "patient_id") |>
    left_join(msss_rr_phase, by = "patient_id")

patients <- patients |>
    left_join(treatment_times, by = "patient_id") |>
    left_join(severity_scores, by = "patient_id")

library(writexl)
dir.create("output", showWarnings = FALSE)
write_xlsx(patients, "output/data.xlsx")
