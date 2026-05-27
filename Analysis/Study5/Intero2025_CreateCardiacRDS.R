# =============================================================================
#  Intero2025_CreateCardiacRDS.R
#
#  Creates per-participant cardiac (heart) RDS files from raw .acq files.
#  Mirrors the structure of the respiration RDS (Intero2025_BeltQualityScreen.R)
#  but targets the heart channel extracted by checkParticipantFiles().
#
#  Channel layout (from checkParticipantFiles / .acq file):
#    R-side recordings: heart = channel 3  (data$channels[[3]])
#    L-side recordings: heart = channel 4  (data$channels[[4]])
#
#  Downsampling:
#    Native SR: 2000 Hz → Target: 500 Hz via safe_decimate(x, q = 4)
#    For q = 4: safe_decimate stages as [4]; actual ratio = 4;
#    effective SR = 2000 / 4 = 500 Hz exactly (no drift correction needed,
#    unlike the breath pipeline where q = 80 stages as [13, 6] → ratio = 78).
#
#  Output:
#    One .rds per participant → Results/heartRDS/<id>.rds
#    QC summary table        → Results/cardiac_rds_qc.csv
#
#  The 14963 / 14962 physio-file ID mismatch from the breath pipeline is
#  handled identically: after the main loop, 14962.rds is copied to 14963.rds
#  if the latter does not already exist.
#
#  Usage:
#    Set FORCE_RECREATE <- TRUE to overwrite existing heartRDS files.
# =============================================================================

rm(list = ls())

# ── Constants ─────────────────────────────────────────────────────────────────
FORCE_RECREATE     <- FALSE   # Set TRUE to overwrite existing heartRDS files

TARGET_CARDIAC_SR  <- 500     # Hz — adequate for R-peak detection (Pan-Tompkins etc.)
NATIVE_SR          <- 2000    # Hz — Biopac acquisition SR

FLATLINE_WIN_S     <- 5       # s  — window length for flatline detection
FLATLINE_SD_THRESH <- 0.01    # SD below this threshold = flatlined window
FLATLINE_FLAG_PCT  <- 20      # % flatlined windows above which we flag the participant

# ── Paths ─────────────────────────────────────────────────────────────────────
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

conditionFile <- paste0(dataPath,    'ConditionLookup.xlsx')
taskDataFile  <- paste0(resultsPath, 'dataFile.csv')

# Source files (provide safe_decimate() and checkParticipantFiles())
pipelineFile  <- file.path(analysisPath, "breath_pipeline.R")
functionsFile <- file.path(analysisPath, "Intero2025_RespirationFunctions.R")

heartRDSDir   <- paste0(resultsPath, 'heartRDS/')
cardiacQCFile <- paste0(resultsPath, 'cardiac_rds_qc.csv')

# ── Set Up ────────────────────────────────────────────────────────────────────
## Load libraries -------
packages <- c(
  "readxl",
  "tidyverse",
  "reticulate",
  "signal"
)
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
options(readr.show_col_types = FALSE)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

reticulate::py_install("bioread")
bioread <- import("bioread")

source(pipelineFile)
source(functionsFile)   # provides checkParticipantFiles()

# Create output directory
if (!dir.exists(heartRDSDir)) dir.create(heartRDSDir, recursive = TRUE)


# ── Helper: compute actual decimation ratio after safe_decimate staging ───────
#  Mirrors compute_actual_q() in Intero2025_BehaviourLedBreathAnalysis.R.
#  For q = 4: stages [4] → actual = 4; effective SR = 2000/4 = 500 Hz.
#  For q = 80: stages [13, 6] → actual = 78; effective SR = 2000/78 ≠ 25 Hz.
#  For the cardiac case (q = 4) the effective SR equals the nominal SR exactly,
#  so no drift correction is needed. Included here for completeness and to
#  guard against future changes to the target SR.
compute_actual_q <- function(q) {
  remaining <- as.integer(q)
  actual    <- 1L
  while (remaining > 1L) {
    step      <- min(remaining, 13L)
    actual    <- actual * step
    remaining <- remaining %/% step
    if (remaining <= 1L) break
  }
  actual
}


# ── Helper: flatline detection ────────────────────────────────────────────────
#  Divides the signal into non-overlapping windows of FLATLINE_WIN_S seconds.
#  Returns the percentage of windows whose SD falls below FLATLINE_SD_THRESH.
#  Returns NA if the signal is shorter than one full window.
detect_flatline_pct <- function(signal, sr,
                                window_s  = FLATLINE_WIN_S,
                                sd_thresh = FLATLINE_SD_THRESH) {
  win_samps <- as.integer(round(sr * window_s))
  n         <- length(signal)
  if (n < win_samps) return(NA_real_)
  n_windows <- floor(n / win_samps)
  sds <- vapply(seq_len(n_windows), function(i) {
    seg <- signal[((i - 1L) * win_samps + 1L):(i * win_samps)]
    sd(seg, na.rm = TRUE)
  }, numeric(1L))
  mean(sds < sd_thresh, na.rm = TRUE) * 100   # % flatlined windows
}


# ── Compute decimation parameters ─────────────────────────────────────────────
q_nominal         <- as.integer(NATIVE_SR / TARGET_CARDIAC_SR)   # 4
q_actual          <- compute_actual_q(q_nominal)                   # 4
ACTUAL_CARDIAC_SR <- NATIVE_SR / q_actual                          # 500.0 Hz

message(sprintf(
  "[Setup] Cardiac decimation: %d Hz → %g Hz  (q_nominal=%d, q_actual=%d, effective SR=%.3f Hz)",
  NATIVE_SR, TARGET_CARDIAC_SR, q_nominal, q_actual, ACTUAL_CARDIAC_SR
))


# ── Load participant list ─────────────────────────────────────────────────────
condData <- read_excel(conditionFile) %>%
  dplyr::distinct() %>%
  dplyr::select(id, Computer, Condition) %>%
  dplyr::rename(firstCondition = Condition)
condData$id <- as.numeric(condData$id)

longData <- read.csv(taskDataFile)
longData <- left_join(longData, condData, by = "id")

longData$currentCondition <-
  factor(longData$firstCondition == "Visual" & longData$ses == 1,
         labels = c("Breath", "Visual"))

# condLookup: required by checkParticipantFiles() to count log files and
# determine first condition.  Must contain columns: id, ses, currentCondition.
condLookup <- longData %>%
  dplyr::group_by(id, ses) %>%
  dplyr::select(currentCondition) %>%
  dplyr::distinct()

# Remove test IDs
testingIds <- c("1234", "5678")
longData   <- longData[!longData$id %in% testingIds, ]

idList   <- unique(longData$id)
fileList <- list.files(path = physioPath)

message(sprintf("[Setup] %d participants to process. heartRDS dir: %s",
                length(idList), heartRDSDir))


# ── Main loop ─────────────────────────────────────────────────────────────────
qcRows <- vector("list", length(idList))

for (thisId in seq_along(idList)) {

  currentId <- idList[thisId]

  message(sprintf(
    "\n========== Cardiac RDS  P#%s  (%d / %d) ==========",
    currentId, thisId, length(idList)
  ))

  rdsPath <- file.path(heartRDSDir, paste0(currentId, ".rds"))

  # Initialise QC row — updated as we learn more about this participant
  qcRow <- data.frame(
    id                = currentId,
    physio_file_found = FALSE,
    side              = NA_character_,
    abcode            = NA_integer_,
    n_triggers        = NA_integer_,
    native_sr         = NATIVE_SR,
    target_sr         = TARGET_CARDIAC_SR,
    actual_sr         = ACTUAL_CARDIAC_SR,
    raw_n_samples     = NA_integer_,
    raw_duration_s    = NA_real_,
    ds_n_samples      = NA_integer_,
    ds_duration_s     = NA_real_,
    signal_mean       = NA_real_,   # computed on raw signal
    signal_sd         = NA_real_,   # computed on raw signal
    signal_range      = NA_real_,   # computed on raw signal
    pct_flatline      = NA_real_,   # % of 5-s windows with SD < FLATLINE_SD_THRESH
    rds_saved         = FALSE,
    rds_path          = rdsPath,
    notes             = "",
    stringsAsFactors  = FALSE
  )

  # ── Skip if already saved (unless FORCE_RECREATE) ─────────────────────────
  if (file.exists(rdsPath) && !FORCE_RECREATE) {
    message(sprintf("  [SKIP] heartRDS already exists: %s", basename(rdsPath)))
    qcRow$rds_saved <- TRUE
    qcRow$notes     <- "Skipped — existing RDS retained"
    qcRows[[thisId]] <- qcRow
    next
  }

  # ── Read .acq file via checkParticipantFiles() ────────────────────────────
  checkResults <- tryCatch(
    checkParticipantFiles(
      currentId, condLookup, fileList, physioPath, bioread,
      expectedLogs         = 2,
      expectedTriggerVal_R = 1,
      expectedTriggerVal_L = 16
    ),
    error = function(e) {
      message(paste("  [ERROR] checkParticipantFiles failed:", e$message))
      NULL
    }
  )

  if (is.null(checkResults)) {
    qcRow$notes <- "checkParticipantFiles crashed"
    qcRows[[thisId]] <- qcRow
    next
  }

  # Populate QC from checkResults (even if heart_data is absent)
  qcRow$side       <- checkResults$side
  qcRow$abcode     <- checkResults$abcode
  qcRow$n_triggers <- checkResults$nTriggers

  # Guard: no heart channel returned
  if (is.null(checkResults$heart_data) || length(checkResults$heart_data) == 0) {
    msg <- sprintf("No heart data returned (abcode=%s, side=%s)",
                   checkResults$abcode, checkResults$side)
    message(sprintf("  [WARN] %s", msg))
    qcRow$notes <- msg
    qcRows[[thisId]] <- qcRow
    next
  }

  qcRow$physio_file_found <- TRUE
  heart_raw               <- checkResults$heart_data
  qcRow$raw_n_samples     <- length(heart_raw)
  qcRow$raw_duration_s    <- length(heart_raw) / NATIVE_SR

  message(sprintf("  Raw heart signal: %d samples @ %d Hz = %.1f s",
                  qcRow$raw_n_samples, NATIVE_SR, qcRow$raw_duration_s))

  # ── Signal quality metrics (computed on raw signal) ───────────────────────
  qcRow$signal_mean  <- mean(heart_raw,  na.rm = TRUE)
  qcRow$signal_sd    <- sd(heart_raw,    na.rm = TRUE)
  qcRow$signal_range <- diff(range(heart_raw, na.rm = TRUE))
  qcRow$pct_flatline <- detect_flatline_pct(heart_raw, NATIVE_SR)

  if (!is.na(qcRow$signal_sd) && qcRow$signal_sd < FLATLINE_SD_THRESH) {
    msg <- "Entire signal flatlined (very low SD)"
    message(sprintf("  [WARN] %s", msg))
    qcRow$notes <- paste(qcRow$notes, "|", msg)
  } else if (!is.na(qcRow$pct_flatline) && qcRow$pct_flatline > FLATLINE_FLAG_PCT) {
    msg <- sprintf("%.1f%% of signal flatlined — check belt continuity", qcRow$pct_flatline)
    message(sprintf("  [WARN] %s", msg))
    qcRow$notes <- paste(qcRow$notes, "|", msg)
  }

  # ── Downsample: 2000 Hz → 500 Hz using safe_decimate(x, q = 4) ───────────
  #  safe_decimate applies an anti-aliasing filter then decimates.
  #  For q = 4 (≤ 13) no multi-stage decomposition is needed, so the
  #  effective SR equals the nominal SR exactly (= 500 Hz).
  heart_ds <- tryCatch(
    safe_decimate(heart_raw, q_nominal),
    error = function(e) {
      message(paste("  [ERROR] safe_decimate failed:", e$message))
      NULL
    }
  )

  if (is.null(heart_ds)) {
    qcRow$notes <- paste(qcRow$notes, "| safe_decimate failed")
    qcRow$notes <- trimws(gsub("^\\| | \\|$", "", qcRow$notes))
    qcRows[[thisId]] <- qcRow
    next
  }

  qcRow$ds_n_samples  <- length(heart_ds)
  qcRow$ds_duration_s <- length(heart_ds) / ACTUAL_CARDIAC_SR

  message(sprintf("  Downsampled:      %d samples @ %.3f Hz = %.1f s",
                  qcRow$ds_n_samples, ACTUAL_CARDIAC_SR, qcRow$ds_duration_s))

  # ── Trigger start indices ─────────────────────────────────────────────────
  #  Stored in NATIVE-SR samples (not downsampled), consistent with the
  #  non-trimmed breath RDS convention.  Load code divides by native_fs to
  #  convert to recording-absolute seconds:
  #    start_times_s = start_indices / native_fs
  start_indices <- checkResults$startIndices   # NULL if no triggers

  if (!is.null(start_indices) && length(start_indices) > 0) {
    message(sprintf("  Trigger indices (native SR): %s",
                    paste(start_indices, collapse = ", ")))
  } else {
    message("  No trigger indices available (abcode = 0 or 1)")
  }

  # ── Build RDS object ──────────────────────────────────────────────────────
  #  Field naming mirrors the breath RDS (from Intero2025_BeltQualityScreen.R)
  #  so that load_participant_physio()-style wrappers can be adapted easily.
  #
  #  Fields:
  #    signal            — downsampled heart signal (500 Hz)
  #    fs                — nominal target SR (500 Hz)
  #    native_fs         — acquisition SR (2000 Hz)
  #    actual_fs         — effective SR after safe_decimate staging (500 Hz)
  #                        For q = 4 this equals fs exactly; included for
  #                        transparency and to match the breath RDS pattern.
  #    start_indices     — trigger positions in NATIVE-SR samples (or NULL)
  #    abcode            — trigger QC code (2=both, 1=one, 0=none, 99/98=error)
  #    n_triggers        — count of valid start triggers found
  #    side              — "L" or "R" (which Biopac port was used)
  #    signal_duration_s — duration of the downsampled signal (s)
  #    raw_duration_s    — duration of the original raw signal (s)
  #    pct_flatline      — % of 5-s windows with SD < 0.01 (QC indicator)
  cardiac_rds <- list(
    signal            = heart_ds,
    fs                = TARGET_CARDIAC_SR,
    native_fs         = NATIVE_SR,
    actual_fs         = ACTUAL_CARDIAC_SR,
    start_indices     = start_indices,
    abcode            = checkResults$abcode,
    n_triggers        = checkResults$nTriggers,
    side              = checkResults$side,
    signal_duration_s = qcRow$ds_duration_s,
    raw_duration_s    = qcRow$raw_duration_s,
    pct_flatline      = qcRow$pct_flatline
  )

  # ── Save RDS ──────────────────────────────────────────────────────────────
  saved_ok <- tryCatch({
    saveRDS(cardiac_rds, rdsPath)
    TRUE
  }, error = function(e) {
    message(paste("  [ERROR] saveRDS failed:", e$message))
    FALSE
  })

  qcRow$rds_saved <- saved_ok

  if (saved_ok) {
    message(sprintf("  [OK] Saved → %s", rdsPath))
  } else {
    qcRow$notes <- paste(qcRow$notes, "| saveRDS failed")
  }

  qcRow$notes  <- trimws(gsub("^\\| | \\|$", "", qcRow$notes))
  qcRows[[thisId]] <- qcRow

}  # end participant loop


# ── ID mismatch fix: 14963 physio file was saved under 14962 ─────────────────
#  Mirrors the same fix in Intero2025_BehaviourLedBreathAnalysis.R.
#  checkParticipantFiles() greps for the participant ID in the filename, so it
#  finds 14962's .acq when processing that ID but finds nothing for 14963.
#  Copying the heartRDS ensures 14963 is treated as having physio data.
hr_14963 <- file.path(heartRDSDir, "14963.rds")
hr_14962 <- file.path(heartRDSDir, "14962.rds")
if (!file.exists(hr_14963) && file.exists(hr_14962)) {
  file.copy(hr_14962, hr_14963)
  message("\n[INFO] Copied 14962.rds → 14963.rds in heartRDS (ID mismatch fix)")
}


# ── Write QC summary ──────────────────────────────────────────────────────────
qcTable <- do.call(rbind, qcRows)
write.csv(qcTable, cardiacQCFile, row.names = FALSE)

message(sprintf("\n[DONE] Cardiac RDS creation complete."))
message(sprintf("  QC summary written to: %s", cardiacQCFile))
message(sprintf("  Participants attempted  : %d", nrow(qcTable)))
message(sprintf("  RDS saved successfully  : %d", sum(qcTable$rds_saved,         na.rm = TRUE)))
message(sprintf("  Physio file not found   : %d", sum(!qcTable$physio_file_found, na.rm = TRUE)))
message(sprintf("  Skipped (existing RDS)  : %d", sum(grepl("Skipped",  qcTable$notes))))
message(sprintf("  Flatline flagged        : %d", sum(grepl("flatline", qcTable$notes))))
message(sprintf("  No triggers (abcode=0)  : %d", sum(!is.na(qcTable$abcode) & qcTable$abcode == 0)))

