library(ggplot2)
library(rlang)
library(stringr)
library(survival)
library(survminer)
library(xfun)

unlink("output", recursive = TRUE)
dir.create("output")

source("src/data.r")

lm_vars <- c(
    "age_at_onset",
    "arr_during_rr_phase",
    "msss_global",
    "armss_global",
    "msss_rr_phase",
    "armss_rr_phase"
)

prop_vars <- c(
    "gender", "disease_course"
)

km_vars <- c(
    "date_of_first_dss_3",
    "date_of_first_dss_6",
    "date_of_iedss_3",
    "date_of_iedss_6",
    "progression_onset",
    "date_of_high_efficacy_trt_start"
)

for (allele in drb1_alleles) {
    allele_s <- if_else(allele < 10, str_c("0", allele), as.character(allele))
    allele_name <- str_c("hla-drb1*", allele_s)
    allele_col <- str_c("drb1_", allele_s, "_count")

    output_dir <- file.path("output", allele_name)
    dir.create(output_dir)

    data_x <- patients |>
        drop_na(all_of(allele_col)) |>
        mutate(
            allele_count = factor(.data[[allele_col]], ordered = TRUE),
            has_allele = allele_count > 0,
            is_homozigous = allele_count == 2
        )

    group_vars <- c(
        "has_allele", "is_homozigous", "allele_count"
    )

    compare_vars <- c(
        "disease_duration",
        "time_on_high_efficacy_treatment",
        "time_on_moderate_efficacy_treatment"
    )

    for (grp in group_vars) {
        grp_values <- unique(data_x[[grp]])
        if (length(grp_values) < 2) next

        sink(file.path(output_dir, "means.txt"), append = TRUE)
        for (y_col in compare_vars) {
            sample_not_valid <- FALSE
            for (gv in grp_values) {
                cat(str_glue("> {y_col} ({grp}_{gv})\n"))
                cat("\n")

                grp_entries <- filter(data_x, .data[[grp]] == gv)
                y_values <- grp_entries[[y_col]]

                if (sum(!is.na(y_values)) < 3) {
                    cat("\n")
                    cat("Too few observations, skipping normality tests...\n")
                    cat("\n")
                    sample_not_valid <- TRUE
                    next
                }

                if (length(unique(y_values)) < 2) {
                    cat("\n")
                    cat("All observations are equal, skipping normality tests...\n")
                    cat("\n")
                    sample_not_valid <- TRUE
                    next
                }

                sw_test <- shapiro.test(y_values)
                print(sw_test)
                cat("\n")

                grp_entries |>
                    ggplot(aes(sample = .data[[y_col]])) +
                    geom_qq() +
                    geom_qq_line()

                ggsave(file.path(
                    output_dir,
                    str_glue("{grp}_{gv}-{y_col}-qq.png")
                ))
            }

            cat(">", y_col, "\n")
            cat("\n")

            if (sample_not_valid) {
                cat("There was an issue with the sample, skipping mean comparison tests...\n")
                cat("\n")
                next
            }

            compare_f <- as.formula(str_glue("{y_col} ~ {grp}"))
            if (length(grp_values) == 2) {
                t.test(compare_f, data = data_x) |>
                    print()
                cat("\n")

                wilcox.test(compare_f, data = data_x) |>
                    print()
                cat("\n")
            } else if (length(grp_values) > 2) {
                aov(compare_f, data = data_x) |>
                    print()
                cat("\n")

                kruskal.test(compare_f, data = data_x) |>
                    print()
                cat("\n")
            }
        }
        sink()

        sink(file.path(output_dir, "props.txt"), append = TRUE)
        for (y_col in prop_vars) {
            xtab <- table(data_x[[grp]], data_x[[y_col]])
            cat(str_glue("> {y_col} (group by {grp})\n"))
            cat("\n")

            chisq.test(xtab) |>
                print()
            cat("\n")

            fisher.test(xtab) |>
                print()
            cat("\n")
        }
        sink()

        sink(file.path(output_dir, "lreg.txt"), append = TRUE)
        for (y_col in lm_vars) {
            cat(str_glue("> {y_col} (regressing on {grp})\n"))
            cat("\n")

            lm_f <- as.formula(str_glue("{y_col} ~ {grp}"))
            if (str_starts(y_col, "arr_")) {
                lm_fit <- MASS::glm.nb(lm_f, data = data_x)
            } else {
                lm_fit <- lm(lm_f, data = data_x)
            }

            lm_fit |>
                summary() |>
                print()
            cat("\n")
        }
        sink()

        sink(file.path(output_dir, "cox.txt"), append = TRUE)
        for (y_col in km_vars) {
            km_f <- as.formula(str_glue("Surv(duration, event) ~ {grp}"))
            km_data <- data_x |>
                mutate(
                    event = !is.na(.data[[y_col]]),
                    duration = if_else(event,
                        (.data[[y_col]] - ms_onset) / dyears(1),
                        (last_clinical_follow_up - ms_onset) / dyears(1)
                    )
                )

            cat(str_glue("> {y_col} (group by {grp})\n"))
            cat("\n\n")

            km_fit <- coxph(km_f, data = km_data)
            km_fit |>
                summary() |>
                print()
            cat("\n")

            ggsurvplot(
                surv_fit(km_f, data = km_data),
                title = y_col, xlab = "Time from onset",
                legend.title = str_to_upper(allele_name),
                pval = TRUE
            )

            ggsave(file.path(
                output_dir,
                str_glue("{grp}-{y_col}-km.png")
            ))
        }
        sink()
    }
}
