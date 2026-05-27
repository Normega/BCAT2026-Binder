# =============================================================================
#  Intero2025_HeartbeatDetection.R
#
#  Computes heartbeat counting accuracy for the 3 HBD (heartbeat detection)
#  intervals embedded at the end of Session 1.
#
#  Pipeline overview
#  -----------------
#  For each participant:
#    1. Load cardiac RDS (500 Hz ECG; created by Intero2025_CreateCardiacRDS.R).
#    2. Determine Session 1 onset in recording-absolute seconds using a
#       3-level priority ladder:
#         Priority 1 — Trigger: start_indices[1] / native_fs from the RDS.
#                       Covers all abcode >= 1 participants, including those
#                       absent from the respiration alignment table.
#         Priority 2 — Alignment recovery: est_onset_ses1_s from the
#                       AlignmentRecovery sheet of qcFile.xlsx.
#                       Covers ~40 abcode=0 participants.
#         Priority 3 — Manual override table.
#    3. Detect R-peaks across the full cardiac signal (Pan-Tompkins inspired):
#         i.   Bandpass 8-20 Hz Butterworth (isolates QRS complex).
#         ii.  Differentiate -> square -> 150 ms moving-average integration.
#         iii. Adaptive threshold (updated every 1 s) + 250 ms refractory period.
#    4. Load the participant's HBD data from hrReport.csv. This file was
#       produced by the PsychoPy extraction step and contains, per participant:
#         * Onset / Offset columns for each of the 3 intervals (25 s, 35 s, 55 s)
#           in PsychoPy within-session seconds. Recording-absolute time is
#           obtained by adding the Session 1 onset:
#             window_start_abs = ses1_onset_s + OnsetXX
#         * Participant-reported beat counts (short25Count, medium35Count,
#           long55Count). The 3 intervals are presented in randomised order;
#           hrReport columns are always keyed by duration, not presentation order.
#    5. Count R-peaks within each HBD window and compute Schandry accuracy:
#         score_i = 1 - |actual_i - reported_i| / actual_i
#       Also report absolute error, proportional error, estimated HR, and
#       a plausibility flag (HR outside 30-200 BPM).
#    6. Write per-interval results and a participant-level summary to
#       Results/heartbeat_detection_results.xlsx.
#
#  Notes on timing
#  ---------------
#  No LAG correction is applied to cardiac windows. The 200 ms belt
#  hardware-latency correction in the respiration pipeline compensates for the
#  mechanical delay of the chest-expansion belt; ECG electrodes have
#  effectively zero signal-transduction delay, and the 25-55 s HBD windows
#  make a 200 ms offset negligible in any case.
#
#  HBD intervals appear at the end of Session 1 only; not acquired in Session 2.
#
#  Schandry (1981) reference:
#    Schandry, R. (1981). Heart Beat Perception and Emotional Experience.
#    Psychophysiology, 18(4), 483-488.
# =============================================================================

rm(list = ls())


# -- Paths --------------------------------------------------------------------
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

heartRDSDir   <- paste0(resultsPath, 'heartRDS/')
qcFile        <- paste0(resultsPath, 'qcFile.xlsx')    # respiration QC (alignment table)
hrReportFile  <- paste0(resultsPath, 'hrReport.csv')   # HBD timing + reported counts
taskDataFile  <- paste0(resultsPath, 'dataFile.csv')   # for participant ID list

hbdOutputFile <- paste0(resultsPath, 'heartbeat_detection_results.xlsx')
hbdQCFile     <- paste0(resultsPath, 'heartbeat_detection_qc.csv')


# -- Constants ----------------------------------------------------------------
CARDIAC_SR  <- 500L    # Hz - downsampled cardiac signal in RDS
NATIVE_SR   <- 2000L   # Hz - original acquisition SR

# HBD interval nominal durations (s).
# NOTE: third interval is 55 s, not 45 s as in original Schandry (1981).
HBD_DURATIONS_S <- c(25, 35, 55)

# Mapping from hrReport column names to interval index (by duration, not order)
HBD_ONSET_COLS  <- c("Onset25",      "Onset35",       "Onset55")
HBD_OFFSET_COLS <- c("Offset25",     "Offset35",      "Offset55")
HBD_COUNT_COLS  <- c("short25Count", "medium35Count", "long55Count")

# R-peak detection (Pan-Tompkins)
RWAVE_BP_LO       <-  8      # Hz - bandpass low cutoff
RWAVE_BP_HI       <- 20      # Hz - bandpass high cutoff
RWAVE_REFRACT_S   <-  0.250  # s  - refractory period (max ~240 BPM)
RWAVE_INTEG_WIN   <-  0.150  # s  - moving-average integration window
RWAVE_THRESH_FRAC <-  0.5    # fraction of max integrated signal for initial threshold
RWAVE_ADAPT_WIN_S <-  1.0    # s  - adaptive threshold update interval

# HR plausibility bounds
MIN_HR_BPM <- 30
MAX_HR_BPM <- 200

# Schandry score floor (clamp to 0 when |error| > actual)
SCHANDRY_FLOOR <- 0


# -- Set Up -------------------------------------------------------------------
packages <- c("readxl", "writexl", "tidyverse", "signal")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
options(readr.show_col_types = FALSE)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}


# -- Manual overrides ---------------------------------------------------------
#  Only ses = 1 entries are used here (HBD is Session 1 only).
#  Mirrors physio_overrides in Intero2025_BehaviourLedBreathAnalysis.R.
hbd_overrides <- data.frame(
  id               = c(13081L, 13081L),
  ses              = c(1L,     2L),
  override_onset_s = c(2211.622, NA_real_),
  exclude_physio   = c(FALSE,   TRUE),
  note             = c(
    "Trigger onset used; ses1 fully precedes belt flatline at t~3300s",
    "Belt flatlined ~14s after ses1 ended; ses2 task window entirely within flatline"
  ),
  stringsAsFactors = FALSE
)


# =============================================================================
#  HELPER FUNCTIONS
# =============================================================================

# -- detect_r_peaks() ---------------------------------------------------------
#  Pan-Tompkins-inspired R-peak detector.
#  Returns a numeric vector of peak times in seconds relative to the start of
#  ecg_signal (recording-absolute when the signal begins at t = 0).
detect_r_peaks <- function(ecg_signal,
                           sr          = CARDIAC_SR,
                           bp_lo       = RWAVE_BP_LO,
                           bp_hi       = RWAVE_BP_HI,
                           refract_s   = RWAVE_REFRACT_S,
                           integ_win_s = RWAVE_INTEG_WIN,
                           thresh_frac = RWAVE_THRESH_FRAC,
                           adapt_win_s = RWAVE_ADAPT_WIN_S) {

  n <- length(ecg_signal)
  if (n < sr * 2) {
    warning("detect_r_peaks: signal shorter than 2 s - returning no peaks")
    return(numeric(0))
  }

  # 1. Bandpass filter (Butterworth, order 4)
  nyq <- sr / 2
  bf  <- tryCatch(
    signal::butter(4, c(bp_lo / nyq, bp_hi / nyq), type = "pass"),
    error = function(e) NULL
  )
  ecg_filt <- if (!is.null(bf)) signal::filtfilt(bf, ecg_signal) else ecg_signal

  # 2. Differentiate -> square
  ecg_sq <- diff(ecg_filt)^2

  # 3. Moving-average integration
  win_samp  <- max(1L, as.integer(round(integ_win_s * sr)))
  kernel    <- rep(1 / win_samp, win_samp)
  ecg_integ <- Re(stats::filter(c(ecg_sq, rep(0, win_samp - 1L)), kernel,
                                method = "convolution", sides = 1))[seq_len(n)]
  ecg_integ[is.na(ecg_integ)] <- 0

  # 4. Adaptive threshold + refractory period peak finding
  refract_samp <- as.integer(round(refract_s   * sr))
  adapt_samp   <- as.integer(round(adapt_win_s * sr))
  threshold    <- thresh_frac * max(ecg_integ, na.rm = TRUE)
  r_samples    <- integer(0)
  last_peak    <- -refract_samp
  i            <- 1L

  while (i <= n) {
    if ((i %% adapt_samp) == 1L) {
      seg_max <- max(ecg_integ[i:min(n, i + adapt_samp - 1L)], na.rm = TRUE)
      if (seg_max > 0) threshold <- max(threshold * 0.9, thresh_frac * seg_max)
    }
    if (ecg_integ[i] > threshold && (i - last_peak) > refract_samp) {
      search_end  <- min(n, i + refract_samp %/% 2L)
      local_max_i <- which.max(ecg_integ[i:search_end]) + i - 1L
      r_samples   <- c(r_samples, local_max_i)
      last_peak   <- local_max_i
      i           <- local_max_i + refract_samp
    } else {
      i <- i + 1L
    }
  }

  (r_samples - 1L) / sr   # 1-based: sample 1 = t = 0 s
}


# -- count_beats_in_window() --------------------------------------------------
count_beats_in_window <- function(r_peak_times_abs, window_start_s, window_end_s) {
  sum(r_peak_times_abs >= window_start_s & r_peak_times_abs <= window_end_s)
}


# -- schandry_score() ---------------------------------------------------------
#  Schandry (1981) accuracy; clamped to [SCHANDRY_FLOOR, 1].
schandry_score <- function(actual, reported) {
  if (is.na(actual) || is.na(reported) || actual == 0) return(NA_real_)
  max(SCHANDRY_FLOOR, 1 - abs(actual - reported) / actual)
}


# =============================================================================
#  LOAD REFERENCE DATA
# =============================================================================

longData   <- read.csv(taskDataFile)
testingIds <- c("1234", "5678")
longData   <- longData[!longData$id %in% testingIds, ]
idList     <- unique(longData$id)

# hrReport: HBD window timing (PsychoPy within-session seconds) + reported counts.
# Onset/Offset columns are keyed by duration (25/35/55), NOT presentation order,
# because the 3 intervals are randomised across participants.
hrReport <- tryCatch(
  read.csv(hrReportFile, stringsAsFactors = FALSE),
  error = function(e) stop("Could not load hrReport.csv: ", e$message)
)
hrReport$id <- as.numeric(hrReport$id)

# Alignment recovery table (ses1 onsets for abcode = 0 participants)
alignTable <- tryCatch(
  read_excel(qcFile, sheet = "AlignmentRecovery"),
  error = function(e) { message("[WARN] AlignmentRecovery sheet unavailable: ", e$message); NULL }
)

message(sprintf("[Setup] %d participants | hrReport: %d rows | alignment table: %s",
                length(idList),
                nrow(hrReport),
                if (!is.null(alignTable)) sprintf("%d rows", nrow(alignTable)) else "unavailable"))


# =============================================================================
#  MAIN LOOP
# =============================================================================

hbd_interval_rows <- list()
hbd_summary_rows  <- list()

for (thisIdx in seq_along(idList)) {

  currentId <- idList[thisIdx]

  message(sprintf("\n========== HBD P#%s (%d/%d) ==========",
                  currentId, thisIdx, length(idList)))

  sumRow <- data.frame(
    id                   = currentId,
    onset_source         = NA_character_,
    ses1_onset_abs_s     = NA_real_,
    cardiac_rds_found    = FALSE,
    hrreport_found       = FALSE,
    n_intervals_analyzed = 0L,
    mean_schandry        = NA_real_,
    mean_abs_error       = NA_real_,
    mean_hr_bpm          = NA_real_,
    qc_notes             = "",
    stringsAsFactors     = FALSE
  )

  add_note <- function(msg) {
    sumRow$qc_notes <<- trimws(paste(sumRow$qc_notes, "|", msg))
    message(sprintf("  [NOTE] %s", msg))
  }

  # -- 1. LOAD CARDIAC RDS ----------------------------------------------------
  rdsPath     <- file.path(heartRDSDir, paste0(currentId, ".rds"))
  cardiac_rds <- NULL

  if (file.exists(rdsPath)) {
    cardiac_rds <- tryCatch(
      readRDS(rdsPath),
      error = function(e) { add_note(paste("RDS load failed:", e$message)); NULL }
    )
  } else {
    add_note("heartRDS not found")
  }
  sumRow$cardiac_rds_found <- !is.null(cardiac_rds)

  # -- 2. DETERMINE SESSION 1 ONSET -------------------------------------------
  #
  #  Priority 1 - trigger  : start_indices[1] / native_fs from RDS
  #  Priority 2 - recovery : est_onset_ses1_s from AlignmentRecovery sheet
  #  Priority 3 - override : manual override table
  #
  ses1_onset_s <- NA_real_
  onset_source <- "unavailable"

  # Priority 1: trigger from RDS
  if (!is.null(cardiac_rds) &&
      !is.null(cardiac_rds$start_indices) &&
      length(cardiac_rds$start_indices) >= 1) {
    native_fs    <- if (!is.null(cardiac_rds$native_fs)) cardiac_rds$native_fs else NATIVE_SR
    ses1_onset_s <- cardiac_rds$start_indices[1] / native_fs
    onset_source <- "trigger_rds"
    message(sprintf("  Onset (trigger):            %.2f s", ses1_onset_s))
  }

  # Priority 2: alignment recovery
  if (is.na(ses1_onset_s) && !is.null(alignTable)) {
    ar <- alignTable[alignTable$id == currentId, ]
    if (nrow(ar) > 0 &&
        !is.na(ar$onset_source_ses1[1]) &&
        ar$onset_source_ses1[1] %in% c("recovery", "trigger", "manual_override") &&
        !is.na(ar$est_onset_ses1_s[1])) {
      ses1_onset_s <- ar$est_onset_ses1_s[1]
      onset_source <- paste0("alignment_", ar$onset_source_ses1[1])
      message(sprintf("  Onset (alignment recovery): %.2f s", ses1_onset_s))
    }
  }

  # Priority 3: manual override
  ov <- hbd_overrides[hbd_overrides$id == currentId & hbd_overrides$ses == 1, ]
  if (nrow(ov) > 0) {
    if (ov$exclude_physio[1]) {
      ses1_onset_s <- NA_real_
      onset_source <- "excluded_manual"
      add_note(paste("Manual exclusion ses1:", ov$note[1]))
    } else if (!is.na(ov$override_onset_s[1])) {
      ses1_onset_s <- ov$override_onset_s[1]
      onset_source <- "manual_override"
      message(sprintf("  Onset (manual override):    %.2f s", ses1_onset_s))
    }
  }

  sumRow$ses1_onset_abs_s <- ses1_onset_s
  sumRow$onset_source     <- onset_source

  if (is.na(ses1_onset_s)) {
    add_note("No ses1 onset available - skipping participant")
    hbd_summary_rows[[thisIdx]] <- sumRow
    next
  }

  # -- 3. R-PEAK DETECTION ----------------------------------------------------
  r_peak_times_abs <- NULL

  if (!is.null(cardiac_rds) && !is.null(cardiac_rds$signal)) {
    sr  <- if (!is.null(cardiac_rds$actual_fs)) cardiac_rds$actual_fs else
           if (!is.null(cardiac_rds$fs))        cardiac_rds$fs        else CARDIAC_SR
    ecg <- cardiac_rds$signal

    peak_times_rel <- tryCatch(
      detect_r_peaks(ecg, sr = sr),
      error = function(e) { add_note(paste("R-peak detection failed:", e$message)); NULL }
    )

    if (!is.null(peak_times_rel)) {
      r_peak_times_abs <- peak_times_rel   # recording-absolute (signal starts at t = 0)
      message(sprintf("  R-peaks detected: %d over %.0f s recording",
                      length(r_peak_times_abs), length(ecg) / sr))
    }
  } else {
    add_note("No cardiac signal - R-peak detection skipped")
  }

  # -- 4. LOAD HBD DATA FROM hrReport -----------------------------------------
  #
  #  Onset/Offset are PsychoPy within-session seconds.
  #  Recording-absolute window:
  #    window_start_abs = ses1_onset_s + OnsetXX
  #    window_stop_abs  = ses1_onset_s + OffsetXX
  #
  hr_row <- hrReport[hrReport$id == currentId, ]
  sumRow$hrreport_found <- nrow(hr_row) > 0

  if (nrow(hr_row) == 0) {
    add_note("Participant absent from hrReport.csv - cannot compute accuracy")
    hbd_summary_rows[[thisIdx]] <- sumRow
    next
  }

  hr <- hr_row[1, ]

  message(sprintf("  hrReport counts: 25s=%s  35s=%s  55s=%s",
                  hr$short25Count, hr$medium35Count, hr$long55Count))

  # -- 5. COMPUTE ACCURACY PER INTERVAL ---------------------------------------
  interval_rows_this_p <- list()
  analyzed <- 0L

  for (ivl in seq_along(HBD_DURATIONS_S)) {

    dur         <- HBD_DURATIONS_S[ivl]
    pp_start_s  <- suppressWarnings(as.numeric(hr[[ HBD_ONSET_COLS[ivl]  ]]))
    pp_stop_s   <- suppressWarnings(as.numeric(hr[[ HBD_OFFSET_COLS[ivl] ]]))
    rep_beats   <- suppressWarnings(as.numeric(hr[[ HBD_COUNT_COLS[ivl]  ]]))

    win_start_abs <- ses1_onset_s + pp_start_s
    win_stop_abs  <- ses1_onset_s + pp_stop_s
    obs_dur_s     <- pp_stop_s - pp_start_s   # always approx dur +/- 0.03 s

    irow <- data.frame(
      id               = currentId,
      interval_num     = ivl,
      expected_dur_s   = dur,
      obs_dur_s        = obs_dur_s,
      psychopy_start_s = pp_start_s,
      psychopy_stop_s  = pp_stop_s,
      window_start_abs = win_start_abs,
      window_stop_abs  = win_stop_abs,
      reported_beats   = rep_beats,
      actual_beats     = NA_integer_,
      hr_bpm           = NA_real_,
      schandry_score   = NA_real_,
      abs_error        = NA_real_,
      prop_error       = NA_real_,
      hr_plausible     = NA,
      onset_source     = onset_source,
      qc_note          = "",
      stringsAsFactors = FALSE
    )

    if (!is.null(r_peak_times_abs) &&
        !is.na(win_start_abs) && !is.na(win_stop_abs)) {

      n_beats           <- count_beats_in_window(r_peak_times_abs,
                                                 win_start_abs, win_stop_abs)
      irow$actual_beats <- n_beats
      irow$hr_bpm       <- (n_beats / obs_dur_s) * 60
      irow$hr_plausible <- irow$hr_bpm >= MIN_HR_BPM & irow$hr_bpm <= MAX_HR_BPM

      if (n_beats == 0) {
        irow$qc_note      <- "Zero beats detected - check cardiac signal in this window"
        irow$hr_plausible <- FALSE
      } else if (!irow$hr_plausible) {
        irow$qc_note <- sprintf("HR %.1f BPM outside plausible range [%d, %d]",
                                irow$hr_bpm, MIN_HR_BPM, MAX_HR_BPM)
      }

      if (!is.na(rep_beats) && n_beats > 0) {
        irow$schandry_score <- schandry_score(n_beats, rep_beats)
        irow$abs_error      <- abs(n_beats - rep_beats)
        irow$prop_error     <- irow$abs_error / n_beats
      } else if (is.na(rep_beats)) {
        irow$qc_note <- trimws(paste(irow$qc_note, "| No participant response recorded"))
      }

      analyzed <- analyzed + 1L

    } else if (is.null(r_peak_times_abs)) {
      irow$qc_note <- "No R-peaks available"
    } else {
      irow$qc_note <- "Missing window timing"
    }

    irow$qc_note <- trimws(gsub("^\\| | \\|$", "", irow$qc_note))
    interval_rows_this_p[[ivl]] <- irow
  }

  sumRow$n_intervals_analyzed <- analyzed

  if (analyzed > 0) {
    sumRow$mean_schandry  <- mean(sapply(interval_rows_this_p, `[[`, "schandry_score"), na.rm = TRUE)
    sumRow$mean_abs_error <- mean(sapply(interval_rows_this_p, `[[`, "abs_error"),      na.rm = TRUE)
    sumRow$mean_hr_bpm    <- mean(sapply(interval_rows_this_p, `[[`, "hr_bpm"),         na.rm = TRUE)
  }

  message(sprintf("  %d/3 intervals analyzed | mean Schandry = %s",
                  analyzed,
                  if (!is.na(sumRow$mean_schandry)) round(sumRow$mean_schandry, 3) else "NA"))

  hbd_interval_rows <- c(hbd_interval_rows, interval_rows_this_p)
  hbd_summary_rows[[thisIdx]] <- sumRow

}  # end participant loop


# =============================================================================
#  COMPILE AND WRITE OUTPUTS
# =============================================================================

hbd_summary   <- do.call(rbind, Filter(Negate(is.null), hbd_summary_rows))
hbd_intervals <- if (length(hbd_interval_rows) > 0)
  do.call(rbind, hbd_interval_rows) else data.frame()

# -- Console summary ----------------------------------------------------------
message("\n=== FINAL SUMMARY ===")
message(sprintf("Participants processed:         %d", nrow(hbd_summary)))
message(sprintf("  With ses1 onset:              %d", sum(!is.na(hbd_summary$ses1_onset_abs_s))))
message(sprintf("  In hrReport:                  %d", sum(hbd_summary$hrreport_found,          na.rm = TRUE)))
message(sprintf("  >= 1 interval analyzed:       %d", sum(hbd_summary$n_intervals_analyzed > 0, na.rm = TRUE)))
message(sprintf("Overall mean Schandry: %.3f (SD = %.3f)",
                mean(hbd_summary$mean_schandry, na.rm = TRUE),
                sd(hbd_summary$mean_schandry,   na.rm = TRUE)))

message("\nOnset source breakdown:")
print(table(hbd_summary$onset_source, useNA = "ifany"))

# -- Per-interval summary -----------------------------------------------------
if (nrow(hbd_intervals) > 0) {
  plaus_summary <- hbd_intervals %>%
    dplyr::group_by(interval_num, expected_dur_s) %>%
    dplyr::summarise(
      n_with_physio     = sum(!is.na(actual_beats)),
      n_hr_plausible    = sum(hr_plausible,  na.rm = TRUE),
      mean_actual_beats = round(mean(actual_beats,   na.rm = TRUE), 1),
      sd_actual_beats   = round(sd(actual_beats,     na.rm = TRUE), 1),
      mean_reported     = round(mean(reported_beats, na.rm = TRUE), 1),
      mean_schandry     = round(mean(schandry_score, na.rm = TRUE), 3),
      mean_abs_error    = round(mean(abs_error,      na.rm = TRUE), 2),
      .groups = "drop"
    )
  message("\nPer-interval summary:")
  print(as.data.frame(plaus_summary))
}

# -- QC breakdown by onset source ---------------------------------------------
qc_by_source <- hbd_summary %>%
  dplyr::group_by(onset_source) %>%
  dplyr::summarise(
    n                = dplyr::n(),
    n_cardiac_found  = sum(cardiac_rds_found,  na.rm = TRUE),
    n_hrreport_found = sum(hrreport_found,      na.rm = TRUE),
    n_analyzed       = sum(n_intervals_analyzed > 0, na.rm = TRUE),
    mean_schandry    = round(mean(mean_schandry,  na.rm = TRUE), 3),
    sd_schandry      = round(sd(mean_schandry,    na.rm = TRUE), 3),
    mean_abs_error   = round(mean(mean_abs_error, na.rm = TRUE), 2),
    .groups = "drop"
  )

# -- Write outputs ------------------------------------------------------------
write_xlsx(
  list(
    "ParticipantSummary" = hbd_summary,
    "IntervalResults"    = hbd_intervals,
    "OnsetSourceQC"      = qc_by_source
  ),
  hbdOutputFile
)
write.csv(hbd_summary, hbdQCFile, row.names = FALSE)

message(sprintf("\n[DONE] Results written to:\n  %s\n  %s", hbdOutputFile, hbdQCFile))
