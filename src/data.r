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
    left_join(
        edmus_diagnosis |>
            select(patient_id, ms_onset, disease_course, progression_onset),
        by = "patient_id"
    ) |>
    filter(
        !wait_and_see,
        disease_course %in% c("RR", "SP-R", "SP-NR")
    ) |>
    select(-wait_and_see) |>
    mutate(
        disease_duration = last_clinical_follow_up - ms_onset,
    ) |>
    inner_join(
        hla_data |>
            drop_na(codi_edm) |>
            select(codi_edm, starts_with("tipatge_hla_drb1_")),
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
    select(patient_id, date, edss_entered) |>
    group_by(patient_id) |>
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

mssev_params <- patients |>
    left_join(
        clinical_iedss |>
            slice_max(date, by = "patient_id") |>
            select(patient_id, date_of_last_edss = date, iedss = iedss),
        by = "patient_id"
    ) |>
    transmute(
        patient_id,
        edss = iedss,
        ageatedss = (date_of_last_edss - date_of_birth) / dyears(1),
        dd = disease_duration / dyears(1)
    )

patients_armss <- ms_sev(mssev_params, type = "global_armss")$data

patients_msss <- ms_sev(mssev_params, type = "global_msss", omss = TRUE)$data |>
    filter(dd >= 1)

severity_scores <- patients_msss |>
    select(patient_id, iedss = edss, msss = oGMSSS) |>
    left_join(
        patients_armss |> select(patient_id, armss = gARMSS),
        by = "patient_id"
    )

patients <- patients |>
    left_join(severity_scores, by = "patient_id")

library(writexl)
dir.create("output", showWarnings = FALSE)
write_xlsx(patients, "output/data.xlsx")
