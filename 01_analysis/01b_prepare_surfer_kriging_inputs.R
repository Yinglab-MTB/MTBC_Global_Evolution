#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
data_root <- Sys.getenv("MTBC_DATA_DIR", unset = "data")
res_root  <- Sys.getenv("MTBC_RESULTS_DIR", unset = "results")

input_file <- if (length(args) >= 1) args[1] else file.path(data_root, "MTBC_sublineage_country_frequencies_levels_2_3.xlsx")
outdir     <- if (length(args) >= 2) args[2] else file.path(res_root, "surfer_kriging_inputs")

if (!file.exists(input_file)) {
  stop("Input workbook not found: ", input_file)
}
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

sheet_names <- excel_sheets(input_file)
keep_cols <- c("latitude", "longitude", "frequency")

message("Found ", length(sheet_names), " sheets in: ", input_file)

for (sheet in sheet_names) {
  df <- read_excel(input_file, sheet = sheet)
  missing_cols <- setdiff(keep_cols, colnames(df))
  if (length(missing_cols) > 0) {
    message("  [Skip] Sheet [", sheet, "] missing columns: ", paste(missing_cols, collapse = ", "))
    next
  }

  df_out <- df %>%
    select(all_of(keep_cols))

  safe_name <- gsub('[\\\\/:*?"<>|]', "_", sheet)
  sheet_dir <- file.path(outdir, safe_name)
  if (!dir.exists(sheet_dir)) {
    dir.create(sheet_dir, recursive = TRUE)
  }

  output_file <- file.path(sheet_dir, paste0(safe_name, ".xlsx"))
  write.xlsx(df_out, output_file, rowNames = FALSE)
}
