# ============================================================
#  Intero2025_BeltQualityScreen.R
#
#  Lightweight pre-pass: reads each participant's .acq file,
#  extracts the respiratory channel, runs assess_belt_quality(),
#  and writes a summary CSV.
#
#  Optionally saves the downsampled (25 Hz) signal for each
#  participant as an .rds file, so the full pipeline can load
#  from .rds instead of re-reading the large .acq file.
#
#  Run this BEFORE the full runQCLoop(). Use the output CSV
#  to inspect the quality distribution and finalise thresholds
#  in assess_belt_quality() before committing to a cutoff.
#
#  Output files (all written to resultsPath):
#    belt_quality_screen.csv  — one row per participant
#    /rds/<id>.rds            — downsampled signal (if save_rds = TRUE)
#
#  Dependencies:
#    breath_pipeline.R               — safe_decimate(), apply_butter()
#    assess_belt_quality.R           — assess_belt_quality()
#                                      (or paste that function into
#                                       breath_pipeline.R instead)
#    Intero2025_RespirationFunctions.R — checkParticipantFiles()
#    reticulate + bioread            — for reading .acq files
# ============================================================

# Set Up ---------
## Load libraries ---------
packages <- c(
  "readxl", "writexl",
  "tidyverse",
  "signal",
  "reticulate"
)
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
options(readr.show_col_types = FALSE)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}


# Paths -----------------------------------------------------------
# Assumes this script is run from the same environment as
# Intero2025_BehaviourLedBreathAnalysis.R, where mainPath etc. are
# already defined. If running standalone, set these manually.

# -- PATH CONFIGURATION -------------------------------------------------
# Set ROOT_DIR to the repository root before running.
# Default (".") assumes the script is run from the repo root directory.
if (!exists("ROOT_DIR")) ROOT_DIR <- "."
mainPath     <- ROOT_DIR
dataPath     <- file.path(ROOT_DIR, "data")
taskPath     <- file.path(dataPath, "Behaviour")
physioPath   <- file.path(dataPath, "Physio")
analysisPath <- file.path(ROOT_DIR, "study5_processing")
resultsPath  <- file.path(ROOT_DIR, "results")
# -----------------------------------------------------------------------

source(file.path(analysisPath, "breath_pipeline.R"))
source(file.path(analysisPath, "Intero2025_RespirationFunctions.R"))

reticulate::py_install("bioread")
bioread <- reticulate::import("bioread")


# Configuration ---------------------------------------------------

# Set TRUE to save a downsampled .rds signal file per participant.
# Recommended: saves re-reading 2000 Hz .acq files on every future
# pipeline run (~80x file size reduction at 25 Hz).
save_rds <- TRUE

# Belt QC thresholds — provisional, to be tuned after inspecting
# the distribution in belt_quality_screen.csv.
# All arguments are passed through to assess_belt_quality().
# Leave as NULL to use the function's built-in defaults.
qc_args <- list(
  trim_s         = 120,    # trim 2 min from each end
  target_fs      = 25,     # downsample to 25 Hz
  iqr_unusable   = 0.05,
  iqr_degraded   = 0.20,
  band_ratio_min = 0.40,
  freq_lo        = 0.08,
  freq_hi        = 0.35
)


# Load participant list -------------------------------------------
conditionFile <- paste0(dataPath, "ConditionLookup.xlsx")
taskDataFile  <- paste0(resultsPath, "dataFile.csv")

condData <- read_excel(conditionFile) %>%
  dplyr::distinct() %>%
  dplyr::select(id, Computer, Condition) %>%
  dplyr::rename(firstCondition = Condition)
condData$id <- as.numeric(condData$id)

longData <- read.csv(taskDataFile)
longData <- left_join(longData, condData, by = "id")
longData$currentCondition <- factor(
  longData$firstCondition == "Visual" & longData$ses == 1,
  labels = c("Breath", "Visual")
)

condLookup <- longData %>%
  group_by(id, ses) %>%
  dplyr::select(currentCondition) %>%
  distinct()

# Apply same exclusions as main analysis
exclude_ids <- c(1234, 5678, 99998, 99999,                        # test entries
                 5488, 6310, 10105, 10591, 10834,                  # pilot participants
                 11890, 12949, 13498, 13879, 15832,
                 6256, 9982, 13711, 9553)                          # QC exclusions

longData  <- longData[!longData$id %in% exclude_ids, ]
idList    <- unique(longData$id)
fileList  <- list.files(path = physioPath)

# Create output directories
rds_dir <- file.path(resultsPath, "rds")
if (save_rds && !dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)


# Belt quality screening loop -------------------------------------
run_belt_screen <- function(idList, condLookup, fileList,
                            physioPath, bioread,
                            qc_args, save_rds, rds_dir) {

  results <- vector("list", length(idList))

  for (i in seq_along(idList)) {
    currentId <- idList[i]
    message(sprintf("\n[%d/%d] ID %s", i, length(idList), currentId))

    # Skeleton row — filled in as we go
    row <- data.frame(
      id               = currentId,
      physio_file      = NA_character_,
      side             = NA_character_,
      abcode           = NA_integer_,
      n_triggers       = NA_integer_,
      duration_s       = NA_real_,
      usable_duration_s = NA_real_,
      quality          = NA_character_,
      filtered_iqr     = NA_real_,
      band_ratio       = NA_real_,
      dominant_freq_hz = NA_real_,
      flatline_pct     = NA_real_,
      rds_saved        = FALSE,
      screen_note      = "",
      stringsAsFactors = FALSE
    )

    # Step 1: file check + channel extraction
    checkResults <- tryCatch(
      checkParticipantFiles(
        currentId, condLookup, fileList, physioPath, bioread,
        expectedLogs         = 2,
        expectedTriggerVal_R = 1,
        expectedTriggerVal_L = 16
      ),
      error = function(e) {
        message(sprintf("  [ERROR] checkParticipantFiles: %s", e$message))
        NULL
      }
    )

    if (is.null(checkResults)) {
      row$screen_note <- "checkParticipantFiles crashed"
      results[[i]]    <- row
      next
    }

    row$physio_file <- if (!is.na(checkResults$fileName)) checkResults$fileName else ""
    row$side        <- if (!is.na(checkResults$side))     checkResults$side     else ""
    row$abcode      <- checkResults$abcode
    row$n_triggers  <- checkResults$nTriggers

    if (is.null(checkResults$breath_data)) {
      row$screen_note <- "No breath data (file missing or unreadable)"
      results[[i]]    <- row
      next
    }

    # Step 2: belt quality assessment
    bq <- tryCatch(
      do.call(assess_belt_quality,
              c(list(signal = checkResults$breath_data,
                     fs     = checkResults$SR),
                qc_args)),
      error = function(e) {
        message(sprintf("  [ERROR] assess_belt_quality: %s", e$message))
        NULL
      }
    )

    if (!is.null(bq)) {
      row$quality           <- bq$quality
      row$filtered_iqr      <- bq$filtered_iqr
      row$band_ratio        <- bq$band_ratio
      row$dominant_freq_hz  <- bq$dominant_freq_hz
      row$flatline_pct      <- bq$flatline_pct
      row$duration_s        <- bq$duration_s
      row$usable_duration_s <- bq$usable_duration_s
      row$screen_note       <- bq$qc_note
    } else {
      row$screen_note <- "assess_belt_quality failed"
    }

    # Step 3: optionally save downsampled signal as .rds
    # The downsampled signal here comes from safe_decimate() inside
    # assess_belt_quality(). We re-downsample from the FULL (untrimmed)
    # recording so the .rds file retains the complete time series
    # (alignment recovery needs the recording from t=0, not trimmed).
    if (save_rds) {
      rds_path <- file.path(rds_dir, paste0(currentId, ".rds"))
      tryCatch({
        nativeSR  <- checkResults$SR
        target_fs <- if (!is.null(qc_args$target_fs)) qc_args$target_fs else 25
        q         <- floor(nativeSR / target_fs)
        ds_signal <- if (q >= 2) safe_decimate(checkResults$breath_data, q) else
          as.numeric(checkResults$breath_data)
        ds_fs     <- nativeSR / q

        saveRDS(list(
          id         = currentId,
          signal     = ds_signal,
          fs         = ds_fs,
          native_fs  = nativeSR,
          side       = checkResults$side,
          abcode     = checkResults$abcode,
          n_triggers = checkResults$nTriggers,
          start_indices = checkResults$startIndices
        ), rds_path)

        row$rds_saved <- TRUE
        message(sprintf("  [RDS] Saved: %s  (%.0fHz, %.0f samples)",
                        basename(rds_path), ds_fs, length(ds_signal)))
      }, error = function(e) {
        message(sprintf("  [WARN] RDS save failed: %s", e$message))
      })
    }

    results[[i]] <- row
  }

  do.call(rbind, results)
}


# Run -------------------------------------------------------------
message(sprintf("\nBelt quality screen: %d participants\n", length(idList)))

beltScreen <- run_belt_screen(
  idList      = idList,
  condLookup  = condLookup,
  fileList    = fileList,
  physioPath  = physioPath,
  bioread     = bioread,
  qc_args     = qc_args,
  save_rds    = save_rds,
  rds_dir     = rds_dir
)


# Save results ----------------------------------------------------
outFile <- file.path(resultsPath, "belt_quality_screen.csv")
write.csv(beltScreen, outFile, row.names = FALSE)
message(sprintf("\nResults written to: %s", outFile))


# Quick summary ---------------------------------------------------
message("\n=== Quality Distribution ===")
print(table(beltScreen$quality, useNA = "ifany"))

message("\n=== IQR by quality label ===")
print(
  beltScreen %>%
    group_by(quality) %>%
    summarise(
      n         = n(),
      iqr_min   = min(filtered_iqr, na.rm = TRUE),
      iqr_med   = median(filtered_iqr, na.rm = TRUE),
      iqr_max   = max(filtered_iqr, na.rm = TRUE),
      br_med    = median(band_ratio, na.rm = TRUE),
      freq_med  = median(dominant_freq_hz, na.rm = TRUE),
      .groups   = "drop"
    )
)

message("\n=== Abcode breakdown ===")
print(table(beltScreen$abcode, useNA = "ifany"))
