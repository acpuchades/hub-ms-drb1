library(readxl)

source("src/common.r")

hla_data <- read_excel("data/PACIENTES-DNA-def.xls") |>
    select(-1) |>
    normalize_names() |>
    rename(apellidos = "nombre_paciente", nombre = "_5") |>
    mutate(across(starts_with("tipatge_hla_drb1_"), as.integer))
