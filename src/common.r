library(dplyr)
library(stringr)

normalize_names <- function(data) {
    data |>
        rename_with(\(x) x |>
            str_to_lower() |>
            str_replace("^([^a-z])", "_\\1") |>
            str_replace_all("(?:[^a-z0-9]|\\s)+", "_") |>
            str_replace("_*$", ""))
}
