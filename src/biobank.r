library(dplyr)
library(readxl)

source("src/common.r")

biobank_load <- function(path) {
    read_excel(path, guess_max = 9999) |>
        biobank_normalize_names()
}

biobank_normalize_names <- function(data) {
    data |>
        normalize_names() |>
        rename_with(\(x) str_replace(x, "^fc_", ""))
}

biobank_fclinica <- biobank_load("data/biobank-fclinica.xlsx") |>
    mutate(across(c(nhc, sap), as.character))

biobank_nhc2sap <- function(nhc) {
    tibble(nhc = nhc) |>
        left_join(biobank_fclinica |> select(nhc, sap), by = "nhc", multiple = "first") |>
        pull(sap)
}

biobank_edmus2sap <- function(id_edmus_local) {
    tibble(n_edmus = id_edmus_local) |>
        left_join(biobank_fclinica |> select(n_edmus, sap), by = "n_edmus", multiple = "first") |>
        pull(sap)
}
