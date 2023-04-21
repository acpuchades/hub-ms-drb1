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
    filter(!wait_and_see) |>
    select(-wait_and_see) |>
    mutate(
        age_at_onset = (ms_onset - date_of_birth) / dyears(1),
        disease_duration = (last_clinical_follow_up - ms_onset) / dyears(1)
    ) |>
    inner_join(
        hla_data |>
            drop_na(codi_edm) |>
            select(codi_edm, starts_with("tipatge_hla_drb1_")) |>
            rename_with(~ str_replace(.x, "^tipatge_hla_drb1_", "drb1_allele_")) |>
            mutate(drb1_haplotype = str_c(
                pmin(drb1_allele_1, drb1_allele_2),
                pmax(drb1_allele_1, drb1_allele_2),
                sep = "/"
            )),
        by = c(edmus_local_id = "codi_edm"), multiple = "first"
    )

dss_reached <- patients |>
    select(patient_id) |>
    cross_join(tibble(dss_reached = 1:10)) |>
    left_join(
        edmus_clinical |>
            select(patient_id, date, dss_entered) |>
            drop_na() |>
            group_by(patient_id) |>
            arrange(date, .by_group = TRUE) |>
            mutate(dss_reached = cummax(dss_entered)) |>
            ungroup() |>
            select(patient_id, dss_reached, date),
        by = c("patient_id", "dss_reached")
    ) |>
    slice_min(date, by = c("patient_id", "dss_reached"), n = 1, with_ties = FALSE) |>
    group_by(patient_id) |>
    arrange(dss_reached, .by_group = TRUE) |>
    fill(date, .direction = "up") |>
    ungroup()

patients <- patients |>
    left_join(
        dss_reached |> pivot_wider(
            names_from = "dss_reached",
            names_prefix = "date_of_first_dss_",
            values_from = "date"
        ),
        by = "patient_id"
    ) |>
    relocate(starts_with("date_of_first_dss_"), .before = "date_of_iedss_1")

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
        duration_of_rr_phase = case_match(
            disease_course,
            "RR" ~ last_clinical_follow_up - ms_onset,
            "SP" ~ progression_onset - ms_onset,
            "PP" ~ dyears(0)
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
    drop_na() |>
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

drb1_alleles <- patients$drb1_allele_1 |>
    union(patients$drb1_allele_2) |>
    na.omit() |>
    sort()

allele_info <- patients |>
    select(patient_id, starts_with("drb1_allele_")) |>
    cross_join(tibble(allele = drb1_alleles)) |>
    transmute(
        patient_id, allele,
        count = (drb1_allele_1 == allele) + (drb1_allele_2 == allele)
    ) |>
    arrange(patient_id, allele) |>
    mutate(
        allele_s = if_else(allele < 10, str_c("0", allele), as.character(allele))
    ) |>
    pivot_wider(
        id_cols = patient_id,
        names_from = "allele_s",
        names_glue = "drb1_{allele_s}_count",
        values_from = "count"
    ) |>
    drop_na()

time_to_second_line_trt <- edmus_trt_dm |>
    semi_join(patients, by = "patient_id") |>
    filter(inn %in% edmus_trt_dm_high_efficacy_inn) |>
    slice_min(onset_date, by = "patient_id", n = 1, with_ties = FALSE) |>
    select(patient_id, date_of_high_efficacy_trt_start = onset_date)

patients <- patients |>
    left_join(treatment_times, by = "patient_id") |>
    left_join(time_to_second_line_trt, by = "patient_id") |>
    left_join(severity_scores, by = "patient_id") |>
    left_join(allele_info, by = "patient_id") |>
    relocate(matches("drb1_[0-9]{2}_count"), .after = "drb1_haplotype")

library(writexl)
dir.create("output", showWarnings = FALSE)
write_xlsx(patients, "output/data.xlsx")
