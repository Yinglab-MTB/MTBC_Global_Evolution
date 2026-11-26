#!/usr/bin/env Rscript
# ==============================================================================
# Step 01b: Export Country-Level Latitude, Longitude & Sublineage Frequencies
#           for Ordinary Kriging Interpolation (Golden Software Surfer)
#
# Usage:
#   Rscript 01_analysis/01b_prepare_surfer_kriging_inputs.R \
#       [input_excel_path] [output_directory]
#
# Defaults:
#   input_excel_path : <MTBC_DATA_DIR or ./data>/Fig. 2/二级谱系分布频率.xlsx
#   output_directory : <MTBC_RESULTS_DIR or ./results>/surfer_kriging_inputs
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
data_root <- Sys.getenv("MTBC_DATA_DIR", unset = "data")
res_root  <- Sys.getenv("MTBC_RESULTS_DIR", unset = "results")

input_file <- if (length(args) >= 1) args[1] else file.path(data_root, "Fig. 2", "二级谱系分布频率.xlsx")
outdir     <- if (length(args) >= 2) args[2] else file.path(res_root, "surfer_kriging_inputs")

if (!file.exists(input_file)) {
  stop("Input workbook not found: ", input_file)
}
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

sheet_names <- excel_sheets(input_file)
keep_cols <- c("latitude...7", "longitude...8", "frequency")

message("Found ", length(sheet_names), " sheets in: ", input_file)

for (sheet in sheet_names) {
  df <- read_excel(input_file, sheet = sheet)
  missing_cols <- setdiff(keep_cols, colnames(df))
  if (length(missing_cols) > 0) {
    message("  [Skip] Sheet [", sheet, "] missing columns: ", paste(missing_cols, collapse = ", "))
    next
  }

  df_out <- df %>%
    select(all_of(keep_cols)) %>%
    rename(
      latitude  = `latitude...7`,
      longitude = `longitude...8`
    )

  safe_name <- gsub('[\\\\/:*?"<>|]', "_", sheet)
  sheet_dir <- file.path(outdir, safe_name)
  if (!dir.exists(sheet_dir)) {
    dir.create(sheet_dir, recursive = TRUE)
  }

  output_file <- file.path(sheet_dir, paste0(safe_name, ".xlsx"))
  write.xlsx(df_out, output_file, rowNames = FALSE)
}

message("All sublineage Kriging input tables exported to: ", outdir)
