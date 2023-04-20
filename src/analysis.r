library(ggplot2)
library(rlang)
library(stringr)
library(survival)
library(survminer)
library(xfun)

source("src/data.r")

lm_vars <- c(
    "arr_during_rr_phase",
    "msss_global", "armss_global",
    "msss_rr_phase", "armss_rr_phase"
)

km_vars <- c(
    "date_of_iedss_3",
    "date_of_iedss_6",
    "progression_onset",
    "date_of_high_efficacy_trt_start"
)

unlink("output", recursive = TRUE)
dir.create("output")

for (allele in hla_drb1_alleles) {
    allele_s <- if_else(allele < 10, str_c("0", allele), as.character(allele))
    allele_name <- str_c("hla-drb1*", allele_s)
    allele_col <- str_c("drb1_", allele_s, "_count")

    output_dir <- file.path("output", allele_name)
    dir.create(output_dir)
    sink(file.path(output_dir, "results.txt"))

    data_x <- patients |>
        drop_na(allele_col) |>
        mutate(
            has_allele = .data[[allele_col]] > 0,
            is_homozigous = .data[[allele_col]] == 2
        )

    group_cols <- c(
        "has_allele", "is_homozigous"
    )

    for (grp in group_cols) {
        for (y_col in lm_vars) {
            cat(">", y_col, "\n")
            lm_formula <- as.formula(str_c(y_col, "~", grp, sep = ""))
            if (str_starts(y_col, "arr_")) {
                lm_fit <- MASS::glm.nb(lm_formula, data = data_x)
            } else {
                lm_fit <- lm(lm_formula, data = data_x)
            }

            lm_fit |>
                summary() |>
                print()
            cat("\n")
        }

        for (y_col in km_vars) {
            km_formula <- as.formula(str_c("Surv(duration, event) ~ ", grp))
            km_data <- data_x |>
                mutate(
                    event = !is.na(.data[[y_col]]),
                    duration = if_else(event,
                        (.data[[y_col]] - ms_onset) / dyears(1),
                        (last_clinical_follow_up - ms_onset) / dyears(1)
                    )
                )

            cat(">", y_col, "\n\n")
            km_fit <- coxph(km_formula, data = km_data)
            km_fit |>
                summary() |>
                print()
            cat("\n")

            ggsurvplot(
                surv_fit(km_formula, data = km_data),
                title = y_col, xlab = "Time from onset",
                legend.title = str_to_upper(allele_name),
                pval = TRUE
            )
            output_name <- str_c(grp, y_col, sep = "-") |> with_ext(".png")
            ggsave(file.path(output_dir, output_name))
        }
    }
    sink()
}
