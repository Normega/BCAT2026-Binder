# ============================================================
# run_all_cleaning.R — Master Cleaning Script
# Interoception 2025 Paper
#
# Runs all five study cleaning scripts in order.
# Each script reads raw data from base_dir/RawData/ and writes
# cleaned outputs to base_dir/Data/
#
# Study order (new numbering):
#   Study 1: Two-task interoception (dissertation Ch. 3)     → study1_*.csv
#   Study 2: BART impulsivity (dissertation Ch. 4)           → study2_*.csv
#   Study 3: Misattribution of arousal (dissertation Ch. 2)  → study3_*.csv
#   Study 4: Online BCAT prep study (2024)           → study4_*.csv
#   Study 5: Fall 2025 in-person BCAT                → study5_*.csv
#
# Raw input files expected in base_dir/RawData/:
#   Study 1: study1_data.xls
#   Study 2: study3_uncleaned.xls, Study3_demographics.xlsx
#   Study 3: study3_taskdata.xlsx, study3_questionnairedata.xlsx
#   Study 4: study4_longData.csv, study4_testData.csv, study4_questionnairedata.xlsx
#   Study 5: study5_longData.csv, study5_testData.csv, study5_questionnairedata_clean.csv
#
# Output files written to base_dir/Data/:
#   studyX_long.csv        trial-level data
#   studyX_summary.csv     one row per participant
#   studyX_exclusions.csv  exclusion log
#   (studyX_test.csv also written for Studies 4 and 5)
# ============================================================

base_dir = "I:/Shared drives/Interoception 2025/Paper/"
#base_dir   <- "." #insert base directory here
script_dir <- file.path(base_dir, "analysis", "DataCleaning")

data_dir <- file.path(base_dir, "Data")
raw_dir  <- file.path(base_dir, "RawData")  # not included; see OSF
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("Interoception 2025 — Data Cleaning Run")
message("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("========================================\n")

scripts <- c(
  "study1_clean.R",
  "study2_clean.R",
  "study3_clean.R",
  "study4_clean.R",
  "study5_clean.R"
)

for (s in scripts) {
  message("\n--- Running ", s, " ---")
  tryCatch(
    source(file.path(script_dir, s), echo = FALSE),
    error = function(e) {
      message("ERROR in ", s, ": ", conditionMessage(e))
      message("Continuing with next study...")
    }
  )
}

message("\n========================================")
message("All cleaning scripts complete.")
message("Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("========================================")

# Quick output file check
data_dir <- file.path(base_dir, "data")
if (dir.exists(data_dir)) {
  files <- list.files(data_dir, pattern = "study[1-5]_.*\\.csv")
  message("\nFiles written to ", data_dir, ":")
  for (f in sort(files)) {
    full <- file.path(data_dir, f)
    sz   <- file.info(full)$size
    message("  ", f, "  (", format(sz, big.mark = ","), " bytes)")
  }
}
