# ============================================================
#  Intero2025_BreathingAdherence.R
#
#  Extracts trial-level breathing direction and adherence
#  metrics from the respiratory belt pipeline, then models
#  adherence as a predictor of staircase and test trial
#  performance.
#
#  CONCEPT
#  -------
#  On each paced breathing trial the pacer demands a specific
#  trajectory: acceleration (Faster), deceleration (Slower),
#  or steady rate (NoChange).  The stretch belt provides an
#  independent record of what the participant actually did.
#
#  From the 4-5 detected troughs per trial we compute:
#
#  ibi_slope   — linear slope of inter-trough intervals (IBIs)
#                across the trial (seconds per breath position).
#                Negative = breaths getting shorter (accel).
#                Positive = breaths getting longer (decel).
#                Requires >= 3 troughs (>= 2 IBIs).
#
#  direction_match — whether ibi_slope sign agrees with the
#                    trial's demanded direction:
#                    Change < 0 → Faster → expect ibi_slope < 0
#                    Change > 0 → Slower → expect ibi_slope > 0
#                    Change = 0 → NoChange → no demand
#                    Coded: 1 = match, -1 = mismatch, 0 = NoChange
#
#  ibi_slope_std — ibi_slope standardised within participant ×
#                  session, for comparability across participants
#                  with different resting rates.
#
#  These are joined to the behavioural data (Accuracy, Arousal,
#  Confidence) and used to predict trial outcomes — primarily
#  as an index of moment-to-moment mind wandering / adherence.
#
#  PIPELINE
#  --------
#  For each participant with usable alignment:
#    1. Load downsampled .rds (25.641 Hz effective)
#    2. Clip artifacts → bandpass → z-score → detect troughs
#       (identical parameters to main QC loop)
#    3. For each staircase and test trial, extract troughs
#       within the trial window using the recovered onset
#    4. Compute IBI sequence and linear slope
#    5. Derive direction_match and continuous adherence score
#
#  OUTPUT
#  ------
#  breathAdherence.csv  — one row per trial (staircase + test),
#                         all participants, with IBI metrics
#                         joined to behavioural outcomes.
#  adherence_summary.csv — participant × session summary.
#
#  MODELS  (see Analysis section at bottom)
#  ------
#  H-physio-A  : direction_match → Accuracy (controls: Change, Change²)
#  H-physio-B  : falsification in Visual ses1 (no pacing)
#  H-physio-C  : ibi_slope_std × Change interaction on Accuracy
#  H-physio-D  : direction_match → Arousal
#  H-physio-E  : cumulative adherence (running mean within session) → Threshold
# ============================================================


# Set Up ---------
## Load libraries ---------
packages <- c(
  "readxl", "writexl",
  "tidyverse",
  "signal",
  "lme4", "lmerTest",
  "ggeffects", "patchwork"
)
# Note: pracma not required — trough detection uses run_pipeline() throughout
# for consistency with the main QC pipeline and alignment recovery.
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
options(readr.show_col_types = FALSE)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}


# Paths -----------------------------------------------------------
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

source(file.path(analysisPath, "breath_pipeline.R"))   # run_pipeline(), compute_actual_q()

rds_dir         <- file.path(resultsPath, "rds")
qcFile          <- file.path(resultsPath, "qcFile.xlsx")
taskDataFile    <- file.path(resultsPath, "dataFile.csv")
taskTestFile    <- file.path(resultsPath, "testFile.csv")
outFile         <- file.path(resultsPath, "breathAdherence.csv")
summaryOutFile  <- file.path(resultsPath, "adherence_summary.csv")


# Constants -------------------------------------------------------
STARTDUR   <- 4      # base breath duration (s)
NUMBREATHS <- 4      # breaths per trial
LAG        <- 0.2    # hardware latency correction (s)
MIN_TROUGHS_FOR_SLOPE <- 3   # minimum troughs needed to compute IBI slope


# ── Helper: IBI slope extraction ─────────────────────────────────────────────
#
# Given a vector of trough times (recording-absolute seconds) within
# a trial window, returns a list with:
#   $ibi          — inter-trough intervals (seconds), length n_troughs - 1
#   $ibi_slope    — linear slope of IBI across breath positions (s / breath)
#   $n_ibis       — number of IBIs available
#   $mean_bpm     — mean breathing rate across trial (bpm)
#   $note         — character flag if computation was skipped
#
# Direction convention (matches Change sign):
#   ibi_slope < 0  →  breathing accelerated  (consistent with Change < 0)
#   ibi_slope > 0  →  breathing decelerated  (consistent with Change > 0)
# ─────────────────────────────────────────────────────────────────────────────
compute_ibi_slope <- function(trough_times_abs,
                              trial_start_s,
                              trial_stop_s,
                              pad_s = 1.0) {

  # Filter to padded window
  in_win <- trough_times_abs >= (trial_start_s - pad_s) &
            trough_times_abs <= (trial_stop_s  + pad_s)
  tr     <- sort(trough_times_abs[in_win])

  if (length(tr) < MIN_TROUGHS_FOR_SLOPE) {
    return(list(ibi       = numeric(0),
                ibi_slope = NA_real_,
                n_ibis    = 0L,
                mean_bpm  = NA_real_,
                note      = sprintf("only %d trough(s) — slope not computed",
                                    length(tr))))
  }

  ibi       <- diff(tr)             # inter-trough intervals (s)
  n_ibis    <- length(ibi)
  mean_bpm  <- 60 / mean(ibi)

  # Linear regression of IBI on breath position (1, 2, …, n_ibis)
  # slope < 0: IBIs shrinking → breathing faster
  # slope > 0: IBIs growing  → breathing slower
  positions <- seq_len(n_ibis)
  ibi_slope <- if (n_ibis >= 2) {
    coef(lm(ibi ~ positions))[["positions"]]
  } else {
    NA_real_   # only 1 IBI: can't fit slope
  }

  list(ibi       = ibi,
       ibi_slope = ibi_slope,
       n_ibis    = n_ibis,
       mean_bpm  = mean_bpm,
       note      = "")
}


# ── Helper: direction match ───────────────────────────────────────────────────
#
# Returns 1 (match), -1 (mismatch), or 0 (NoChange trial or insufficient data).
# Match means: sign(ibi_slope) agrees with the direction demanded by Change.
#   Change < 0 → Faster → expect ibi_slope < 0
#   Change > 0 → Slower → expect ibi_slope > 0
# ─────────────────────────────────────────────────────────────────────────────
compute_direction_match <- function(ibi_slope, change) {
  if (is.na(ibi_slope) || change == 0) return(0L)
  if (sign(ibi_slope) == sign(change)) 1L else -1L
}


# Preprocessing is handled by run_pipeline() — no separate helper needed.
# run_pipeline() clips artifacts, bandpass filters, z-scores, and detects
# troughs using the same prominence-based algorithm as the main QC pipeline,
# ensuring trough times are identical across the QC and adherence extractions.


# ── Load data ─────────────────────────────────────────────────────────────────
message("Loading behavioural and QC data...")

longData <- read.csv(taskDataFile)
longTest  <- read.csv(taskTestFile)

# Apply standard exclusions
exclude_ids <- c(1234, 5678, 99998, 99999,
                 5488, 6310, 10105, 10591, 10834,
                 11890, 12949, 13498, 13879, 15832,
                 6256, 9982, 13711, 9553)
longData <- longData[!longData$id %in% exclude_ids, ]
longTest  <- longTest[ !longTest$id  %in% exclude_ids, ]

# QC file: alignment table for onset times and belt quality
alignTable <- read_excel(qcFile, sheet = "AlignmentRecovery")
qcSummary  <- read_excel(qcFile, sheet = "ResultsSummary")

# Normalise test data to staircase column format (matches alignSignals helper)
# Salience: 0 → "Low", 1 → "High"
# Change: Direction × level (signed)
longTest$Salience_str <- ifelse(longTest$Salience == 1, "High", "Low")
longTest$Change       <- longTest$Direction * longTest$level

# Sequence trial numbers within participant × session
longData <- longData %>%
  dplyr::group_by(id, ses) %>%
  dplyr::mutate(trial_seq = dplyr::row_number()) %>%
  dplyr::ungroup()

idList <- unique(longData$id)


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN EXTRACTION LOOP
# ══════════════════════════════════════════════════════════════════════════════
message(sprintf("Extracting breathing adherence for %d participants...", length(idList)))

all_rows <- vector("list", length(idList))

for (p_idx in seq_along(idList)) {

  pid <- idList[p_idx]
  message(sprintf("\n[%d/%d] P%s", p_idx, length(idList), pid))

  # ── Load alignment info ──────────────────────────────────────────────────
  align_row <- alignTable[alignTable$id == pid, ]
  if (nrow(align_row) == 0) {
    message("  No alignment row — skipping")
    next
  }

  first_condition <- align_row$first_condition[1]
  belt_quality    <- align_row$belt_quality[1]

  # Skip unusable belt
  if (!is.na(belt_quality) && belt_quality == "unusable") {
    message("  Belt unusable — skipping")
    next
  }

  # Onset times: prefer trigger, fall back to recovery estimate
  get_onset <- function(source_col, est_col, trig_col) {
    src <- align_row[[source_col]][1]
    if (!is.na(src) && src == "trigger")   return(align_row[[trig_col]][1])
    if (!is.na(src) && src == "recovery")  return(align_row[[est_col]][1])
    NA_real_
  }

  onset1 <- get_onset("onset_source_ses1", "est_onset_ses1_s", "trigger_onset_ses1_s")
  onset2 <- get_onset("onset_source_ses2", "est_onset_ses2_s", "trigger_onset_ses2_s")

  # ── Load and preprocess RDS ──────────────────────────────────────────────
  rds_path <- file.path(rds_dir, paste0(pid, ".rds"))
  if (!file.exists(rds_path)) {
    message("  No RDS — skipping")
    next
  }

  rds <- tryCatch(readRDS(rds_path),
                  error = function(e) { message("  RDS read failed"); NULL })
  if (is.null(rds) || is.null(rds$signal)) { next }

  # Corrected sampling rate (safe_decimate staging fix)
  native_fs    <- if (!is.null(rds$native_fs)) rds$native_fs else 2000
  q_nominal    <- floor(native_fs / 25)
  q_actual     <- compute_actual_q(q_nominal)
  corrected_fs <- native_fs / q_actual   # 25.641 Hz

  # Run the identical pipeline as the main QC loop.
  # Artifact clipping is done inside run_pipeline via the signal_clipped
  # step, matching the main loop exactly.
  # Since the RDS is already at corrected_fs (~25.641 Hz), the
  # downsample step inside run_pipeline auto-skips (q < 2).
  signal_clipped <- pmax(
    pmin(rds$signal, quantile(rds$signal, 0.995, na.rm = TRUE)),
    quantile(rds$signal, 0.005, na.rm = TRUE))

  pipelineResult <- tryCatch(
    run_pipeline(signal_clipped, corrected_fs,
                 detection_method = "prominence",
                 min_prominence   = 0.15,
                 min_dist_s       = 1,
                 detect_hp_on     = FALSE),
    error = function(e) {
      message(sprintf("  run_pipeline failed: %s", e$message))
      NULL
    }
  )

  if (is.null(pipelineResult) || length(pipelineResult$trough_times) == 0) {
    message("  No troughs detected — skipping")
    next
  }
  trough_times <- pipelineResult$trough_times

  # ── Extract per-trial IBI metrics ────────────────────────────────────────
  pData <- longData[longData$id == pid, ]
  pTest <- longTest[longTest$id == pid, ]

  trial_rows <- vector("list",
                       nrow(pData) + nrow(pTest[!is.na(pTest$TestTrial.started), ]))
  row_idx <- 0L

  for (run in 1:2) {

    run_data  <- pData[pData$ses == run, ]
    run_test  <- pTest[pTest$ses == run & !is.na(pTest$TestTrial.started), ]
    run_onset <- if (run == 1) onset1 else onset2

    # For Visual ses1, record trials but note no pacing expected
    is_paced <- !(first_condition == "Visual" && run == 1)

    if (is.na(run_onset) && is_paced) {
      message(sprintf("  Run %d: no onset — trials logged as NA", run))
    }

    # ── Staircase trials ────────────────────────────────────────────────
    for (i in seq_len(nrow(run_data))) {
      row <- run_data[i, ]
      row_idx <- row_idx + 1L

      base_row <- data.frame(
        id              = pid,
        ses             = run,
        trial_type      = "staircase",
        trial_seq       = row$trial_seq,
        trial_num       = row$Trial,
        condition       = row$Condition,
        salience        = row$Salience,
        direction       = row$Direction,
        change          = row$Change,
        level           = row$level,
        accuracy        = row$Accuracy,
        arousal         = row$Arousal,
        confidence      = row$Confidence,
        first_condition = first_condition,
        belt_quality    = belt_quality,
        is_paced        = is_paced,
        onset_available = !is.na(run_onset),
        # Physio metrics (populated below)
        n_troughs       = NA_integer_,
        n_ibis          = NA_integer_,
        ibi_slope       = NA_real_,
        mean_bpm        = NA_real_,
        direction_match = NA_integer_,
        ibi_1 = NA_real_, ibi_2 = NA_real_,
        ibi_3 = NA_real_, ibi_4 = NA_real_,
        physio_note     = "",
        stringsAsFactors = FALSE
      )

      if (!is.na(run_onset)) {
        start_s <- run_onset + LAG + row$trial.started
        stop_s  <- run_onset + LAG + row$trial.stopped

        ibi_res <- compute_ibi_slope(trough_times, start_s, stop_s)

        base_row$n_troughs       <- length(ibi_res$ibi) + 1L
        base_row$n_ibis          <- ibi_res$n_ibis
        base_row$ibi_slope       <- ibi_res$ibi_slope
        base_row$mean_bpm        <- ibi_res$mean_bpm
        base_row$direction_match <- compute_direction_match(ibi_res$ibi_slope,
                                                            row$Change)
        base_row$physio_note     <- ibi_res$note

        # Store individual IBIs (up to 4 for 5 troughs)
        for (k in seq_len(min(4L, ibi_res$n_ibis)))
          base_row[[paste0("ibi_", k)]] <- ibi_res$ibi[k]
      }

      trial_rows[[row_idx]] <- base_row
    }

    # ── Test trials ─────────────────────────────────────────────────────
    if (nrow(run_test) > 0) {
      for (i in seq_len(nrow(run_test))) {
        row <- run_test[i, ]
        row_idx <- row_idx + 1L

        base_row <- data.frame(
          id              = pid,
          ses             = run,
          trial_type      = "test",
          trial_seq       = NA_integer_,
          trial_num       = row$testTrials.thisN + 1L,
          condition       = row$Condition,
          salience        = row$Salience,
          direction       = row$Direction,
          change          = row$Change,
          level           = row$level,
          accuracy        = row$Accuracy,
          arousal         = row$Arousal,
          confidence      = row$Confidence,
          first_condition = first_condition,
          belt_quality    = belt_quality,
          is_paced        = is_paced,
          onset_available = !is.na(run_onset),
          n_troughs       = NA_integer_,
          n_ibis          = NA_integer_,
          ibi_slope       = NA_real_,
          mean_bpm        = NA_real_,
          direction_match = NA_integer_,
          ibi_1 = NA_real_, ibi_2 = NA_real_,
          ibi_3 = NA_real_, ibi_4 = NA_real_,
          physio_note     = "",
          stringsAsFactors = FALSE
        )

        if (!is.na(run_onset)) {
          start_s <- run_onset + LAG + row$TestTrial.started
          stop_s  <- run_onset + LAG + row$TestTrial.stopped

          ibi_res <- compute_ibi_slope(trough_times, start_s, stop_s)

          base_row$n_troughs       <- length(ibi_res$ibi) + 1L
          base_row$n_ibis          <- ibi_res$n_ibis
          base_row$ibi_slope       <- ibi_res$ibi_slope
          base_row$mean_bpm        <- ibi_res$mean_bpm
          base_row$direction_match <- compute_direction_match(ibi_res$ibi_slope,
                                                              row$Change)
          base_row$physio_note     <- ibi_res$note

          for (k in seq_len(min(4L, ibi_res$n_ibis)))
            base_row[[paste0("ibi_", k)]] <- ibi_res$ibi[k]
        }

        trial_rows[[row_idx]] <- base_row
      }
    }
  } # end run loop

  all_rows[[p_idx]] <- do.call(rbind, Filter(Negate(is.null), trial_rows))

} # end participant loop


# ══════════════════════════════════════════════════════════════════════════════
#  POST-PROCESSING
# ══════════════════════════════════════════════════════════════════════════════
message("\nAssembling output...")

breathData <- do.call(rbind, Filter(Negate(is.null), all_rows))

# ── Within-participant standardisation of ibi_slope ─────────────────────────
# Standardise ibi_slope within participant × session so that individual
# differences in resting breath rate don't confound the direction effect.
breathData <- breathData %>%
  dplyr::group_by(id, ses) %>%
  dplyr::mutate(
    ibi_slope_std = if (sum(!is.na(ibi_slope)) >= 3)
      (ibi_slope - mean(ibi_slope, na.rm = TRUE)) / sd(ibi_slope, na.rm = TRUE)
    else
      NA_real_,
    # Cumulative adherence: rolling mean of direction_match up to this trial
    # (within paced staircase only — excludes NoChange trials and test trials)
    adherence_cumulative = {
      dm <- ifelse(trial_type == "staircase" & change != 0 & is_paced,
                   direction_match, NA_real_)
      # Lagged: use adherence from PRIOR trials to predict current outcome
      dplyr::lag(cummean(ifelse(is.na(dm), 0, dm)))
    }
  ) %>%
  dplyr::ungroup()

# Factor coding
breathData$direction_match_f <- factor(breathData$direction_match,
                                        levels = c(-1, 0, 1),
                                        labels = c("mismatch", "neutral", "match"))


# ── Participant-level summary ────────────────────────────────────────────────
adhereSummary <- breathData %>%
  dplyr::filter(is_paced, trial_type == "staircase", change != 0,
         !is.na(direction_match)) %>%
  dplyr::group_by(id, ses, first_condition, belt_quality) %>%
  dplyr::summarise(
    n_trials_with_slope = sum(!is.na(ibi_slope)),
    pct_match           = mean(direction_match == 1,  na.rm = TRUE) * 100,
    pct_mismatch        = mean(direction_match == -1, na.rm = TRUE) * 100,
    mean_ibi_slope      = mean(ibi_slope,             na.rm = TRUE),
    mean_bpm            = mean(mean_bpm,               na.rm = TRUE),
    .groups = "drop"
  )

message(sprintf("Adherence summary: %d participant-sessions", nrow(adhereSummary)))
message(sprintf("  Overall direction match: %.1f%%",
                mean(adhereSummary$pct_match, na.rm = TRUE)))


# ── Save ────────────────────────────────────────────────────────────────────
write.csv(breathData,    outFile,        row.names = FALSE)
write.csv(adhereSummary, summaryOutFile, row.names = FALSE)
message(sprintf("Saved:\n  %s\n  %s", outFile, summaryOutFile))


# ══════════════════════════════════════════════════════════════════════════════
#  ANALYSIS MODELS
# ══════════════════════════════════════════════════════════════════════════════
# Run after extraction is complete and data inspected.
# Restrict to paced staircase trials with physio data available.
# ══════════════════════════════════════════════════════════════════════════════

# ── Subset: paced staircase, non-NoChange, physio available ──────────────────
stair_paced <- breathData %>%
  dplyr::filter(
    trial_type      == "staircase",
    is_paced,
    change          != 0,
    !is.na(ibi_slope),
    onset_available
  ) %>%
  dplyr::mutate(
    change_sq = change^2,
    # Centre ibi_slope_std for interaction interpretability
    ibi_slope_std_c = scale(ibi_slope_std, center = TRUE, scale = FALSE)[, 1]
  )

# Visual ses1 falsification subset
visual_ses1 <- breathData %>%
  dplyr::filter(ses == 1, first_condition == "Visual",
         trial_type == "staircase", change != 0,
         !is.na(ibi_slope), onset_available) %>%
  dplyr::mutate(change_sq = change^2)

# Test trials: paced only
test_paced <- breathData %>%
  dplyr::filter(trial_type == "test", is_paced,
         !is.na(ibi_slope), onset_available) %>%
  dplyr::mutate(change_sq = change^2)

# ── H-physio-A: direction_match → Accuracy (staircase) ──────────────────────
# Does breathing in the demanded direction predict detection accuracy,
# over and above Change magnitude?
# direction_match coded as numeric (-1 / 0 / 1): treat as continuous
# (or run as ordered factor if preferred).

model_A <- glmer(
  accuracy ~ direction_match + change + change_sq + (1 | id),
  data   = stair_paced,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
summary(model_A)

# ── H-physio-B: falsification (three-tier) ───────────────────────────────────
#
# Three falsification tests in increasing order of inferential strength:
#
# B1 — Visual ses1 (between-person): different participants, no pacing.
#      Weakest: null could reflect sample composition rather than instruction.
#
# B2 — Visual-first ses1 vs ses2 (within-person): SAME participants did
#      unpaced ses1 then paced ses2. Controls for individual differences in
#      breathing patterns, task engagement, and change magnitude exposure.
#      If IBI slope predicts accuracy only in paced ses2 and not unpaced ses1,
#      that directly implicates the breathing instruction.
#
# B3 — NoChange trials within paced sessions: change = 0 so no directional
#      demand. IBI slope should be unrelated to accuracy here, confirming
#      the effect in Model C is direction-specific rather than a general
#      breathing-variability effect.

# B1: Visual ses1 between-person
model_B1 <- glmer(
  accuracy ~ direction_match + change + change_sq + (1 | id),
  data   = visual_ses1,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
cat("=== H-physio-B1: Visual ses1 (between-person) ===\n")
print(summary(model_B1)$coefficients)

# B2: Within-person — Visual-first group, ses1 (unpaced) vs ses2 (paced)
vf_ses1 <- breathData %>%
  dplyr::filter(first_condition == "Visual", ses == 1,
                trial_type == "staircase", change != 0, !is.na(ibi_slope)) %>%
  dplyr::mutate(change_sq = change^2,
                ibi_slope_std_c = scale(ibi_slope_std, center = TRUE,
                                        scale = FALSE)[, 1])

vf_ses2 <- breathData %>%
  dplyr::filter(first_condition == "Visual", ses == 2,
                trial_type == "staircase", change != 0, !is.na(ibi_slope)) %>%
  dplyr::mutate(change_sq = change^2,
                ibi_slope_std_c = scale(ibi_slope_std, center = TRUE,
                                        scale = FALSE)[, 1])

cat(sprintf("\n=== H-physio-B2: Visual-first within-person ===\n"))
cat(sprintf("  ses1 unpaced: %d trials, %d participants\n",
            nrow(vf_ses1), dplyr::n_distinct(vf_ses1$id)))
cat(sprintf("  ses2 paced:   %d trials, %d participants\n",
            nrow(vf_ses2), dplyr::n_distinct(vf_ses2$id)))

# B2a: unpaced ses1
# NOTE: vf_ses1 is the same dataset as visual_ses1 used in B1 — both filter
# to Visual-first participants in ses1. B2a therefore adds no new information
# over B1. The genuine within-person contribution comes from B2b (same
# participants in their paced ses2) and B2c (continuous interaction comparison
# across sessions). model_B2a reuses the fitted B1 object to avoid redundant
# computation.
model_B2a <- model_B1
cat("\nB2a — Visual-first ses1 (unpaced) [identical to B1 — see note]:\n")
print(summary(model_B2a)$coefficients)
cat("  --> Within-person falsification: compare B2b (paced ses2) and B2c\n")

# B2b: paced ses2 — should replicate Model A
model_B2b <- glmer(
  accuracy ~ direction_match + change + change_sq + (1 | id),
  data = vf_ses2, family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
cat("\nB2b — Visual-first ses2 (paced):\n")
print(summary(model_B2b)$coefficients)

# B2c: continuous ibi_slope interaction within Visual-first (replicate Model C)
model_B2c_unpaced <- glmer(
  accuracy ~ ibi_slope_std_c * change + change_sq + (1 | id),
  data = vf_ses1, family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
model_B2c_paced <- glmer(
  accuracy ~ ibi_slope_std_c * change + change_sq + (1 | id),
  data = vf_ses2, family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
cat(sprintf("\nB2c — ibi_slope × change interaction (Visual-first):\n"))
cat(sprintf("  ses1 unpaced β = %.4f  p = %.4f\n",
            fixef(model_B2c_unpaced)["ibi_slope_std_c:change"],
            summary(model_B2c_unpaced)$coefficients[
              "ibi_slope_std_c:change", "Pr(>|z|)"]))
cat(sprintf("  ses2 paced   β = %.4f  p = %.4f\n",
            fixef(model_B2c_paced)["ibi_slope_std_c:change"],
            summary(model_B2c_paced)$coefficients[
              "ibi_slope_std_c:change", "Pr(>|z|)"]))

# B3: NoChange trials within paced sessions
# change = 0: no directional demand; ibi_slope_std should not predict accuracy
nochange_paced <- breathData %>%
  dplyr::filter(trial_type == "staircase", is_paced,
                change == 0, !is.na(ibi_slope)) %>%
  dplyr::mutate(ibi_slope_std_c = scale(ibi_slope_std, center = TRUE,
                                        scale = FALSE)[, 1])

cat(sprintf("\n=== H-physio-B3: NoChange trials within paced sessions ===\n"))
cat(sprintf("  N = %d trials, %d participants\n",
            nrow(nochange_paced), dplyr::n_distinct(nochange_paced$id)))

model_B3 <- glmer(
  accuracy ~ ibi_slope_std_c + (1 | id),
  data = nochange_paced, family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
cat("\nB3 — ibi_slope on NoChange accuracy:\n")
print(summary(model_B3)$coefficients)

# Summary
cat("\n=== FALSIFICATION SUMMARY ===\n")
cat(sprintf("Model A   (paced)               direction_match β = %.4f  p = %.6f\n",
            fixef(model_A)["direction_match"],
            summary(model_A)$coefficients["direction_match","Pr(>|z|)"]))
cat(sprintf("Model B1  (Visual ses1 between) direction_match β = %.4f  p = %.4f\n",
            fixef(model_B1)["direction_match"],
            summary(model_B1)$coefficients["direction_match","Pr(>|z|)"]))
cat(sprintf("Model B2a (Vis-first ses1 unpaced) β = %.4f  p = %.4f\n",
            fixef(model_B2a)["direction_match"],
            summary(model_B2a)$coefficients["direction_match","Pr(>|z|)"]))
cat(sprintf("Model B2b (Vis-first ses2 paced)   β = %.4f  p = %.4f\n",
            fixef(model_B2b)["direction_match"],
            summary(model_B2b)$coefficients["direction_match","Pr(>|z|)"]))

# Keep model_B pointing to B1 for backward compatibility with sensitivity code
model_B <- model_B1

# ── H-physio-C: ibi_slope × Change interaction on Accuracy ──────────────────
# Does the DEGREE of breathing rate change (ibi_slope_std) interact
# with Change magnitude to predict accuracy?
# Prediction: when participants breathe in a more extreme direction,
# larger Change × ibi_slope interaction reflects physiological amplification
# of the interoceptive signal.

model_C <- glmer(
  accuracy ~ ibi_slope_std_c * change + change_sq +
             (1 | id),
  data   = stair_paced,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
summary(model_C)

# ── H-physio-D: direction_match → Arousal ───────────────────────────────────
# Trials on which participants followed the pacer's direction should
# produce stronger arousal responses, because the physiological change
# is genuine rather than absent.

model_D <- lmer(
  arousal ~ direction_match + change + change_sq + (1 | id),
  data = stair_paced
)
summary(model_D)

# ── H-physio-E: cumulative adherence → Staircase threshold ──────────────────
# A participant who follows the pacer better (higher cumulative direction
# match across the session) should have lower detection thresholds,
# because their interoceptive signals are more reliable.
# Uses the LAGGED cumulative adherence (prior trials) to avoid circularity.

# Compute participant-level mean adherence per session
adherence_pp <- breathData %>%
  dplyr::filter(trial_type == "staircase", is_paced,
         change != 0, !is.na(direction_match)) %>%
  dplyr::group_by(id, ses) %>%
  dplyr::summarise(mean_adherence = mean(direction_match == 1, na.rm = TRUE),
            .groups = "drop")

# Join to staircase threshold data
stair_thresh <- longData %>%
  dplyr::filter(Trial >= 8, Direction != "NoChange") %>%
  dplyr::group_by(id, ses) %>%
  dplyr::summarise(mean_threshold = mean(level, na.rm = TRUE),
            mean_accuracy  = mean(Accuracy, na.rm = TRUE),
            .groups = "drop") %>%
  dplyr::left_join(adherence_pp, by = c("id", "ses")) %>%
  dplyr::filter(!is.na(mean_adherence))

model_E <- lmer(
  mean_threshold ~ mean_adherence + (1 | id),
  data = stair_thresh
)
summary(model_E)

# Also for test trial accuracy
test_acc <- breathData %>%
  dplyr::filter(trial_type == "test", is_paced, !is.na(direction_match)) %>%
  dplyr::group_by(id, ses) %>%
  dplyr::summarise(test_accuracy  = mean(accuracy, na.rm = TRUE),
            test_adherence = mean(direction_match == 1, na.rm = TRUE),
            .groups = "drop")

model_E2 <- lmer(
  test_accuracy ~ test_adherence + (1 | id),
  data = test_acc
)
summary(model_E2)


# ── Quick descriptive plots ──────────────────────────────────────────────────
# 1. Direction match rate by trial type and condition group
p1 <- adhereSummary %>%
  ggplot(aes(x = factor(ses), y = pct_match, fill = first_condition)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  scale_fill_manual(values = c("Breath" = "#2c7be5", "Visual" = "#e63946")) +
  labs(title = "Direction match rate by session and group",
       x = "Session", y = "% trials where breathing direction matched demand",
       fill = "Group") +
  theme_minimal(base_size = 11)

# 2. ibi_slope distribution by demanded direction
p2 <- stair_paced %>%
  dplyr::mutate(demanded = ifelse(change < 0, "Faster (accel)", "Slower (decel)")) %>%
  ggplot(aes(x = ibi_slope_std, fill = demanded)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("Faster (accel)" = "#2dc653",
                                "Slower (decel)" = "#f4a261")) +
  facet_wrap(~first_condition + ses,
             labeller = label_both) +
  labs(title = "Observed IBI slope distribution by demanded direction",
       x = "IBI slope (standardised within participant × session)",
       y = "Density", fill = "Demanded direction") +
  theme_minimal(base_size = 10)

# 3. Direction match × accuracy
p3 <- stair_paced %>%
  dplyr::group_by(id, ses, direction_match_f, first_condition) %>%
  dplyr::summarise(acc = mean(accuracy, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = direction_match_f, y = acc, colour = first_condition)) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.8) +
  scale_colour_manual(values = c("Breath" = "#2c7be5", "Visual" = "#e63946")) +
  facet_wrap(~ses, labeller = label_both) +
  labs(title = "Accuracy by breathing direction match",
       x = "Direction match", y = "Mean accuracy", colour = "Group") +
  theme_minimal(base_size = 11)

combined_plots <- p1 / (p2 | p3)
ggsave(file.path(resultsPath, "adherence_plots.png"),
       combined_plots, width = 14, height = 12, dpi = 130)
message("Plots saved.")

# ══════════════════════════════════════════════════════════════════════════════
#  SENSITIVITY ANALYSES
#  Pre-specified prior to examining results (see methods markdown 2026-04-13).
#  Address the potential confound between degraded belt quality, low adherence,
#  and accuracy — particularly in the Visual-first Session 2 cluster.
# ══════════════════════════════════════════════════════════════════════════════

# ── Helper: refit primary models on a restricted subset and compare ───────────
compare_sensitivity <- function(primary_model, sens_data,
                                formula, family = NULL,
                                label = "") {
  sens_model <- if (is.null(family)) {
    lmer(formula, data = sens_data)
  } else {
    glmer(formula, data = sens_data, family = family,
          control = glmerControl(optimizer = "bobyqa"))
  }
  primary_coef <- fixef(primary_model)["direction_match"]
  sens_coef    <- fixef(sens_model)["direction_match"]
  cat(sprintf("\n%s\n  Primary β = %.4f  |  Sensitivity β = %.4f  |  N = %d\n",
              label, primary_coef, sens_coef, nrow(sens_data)))
  print(summary(sens_model)$coefficients)
  invisible(sens_model)
}

# ── Sensitivity 1: Good belt only ────────────────────────────────────────────
cat("\n\n══════════════════════════════════════════════\n")
cat("SENSITIVITY 1 — Good belt quality only\n")
cat("══════════════════════════════════════════════\n")

stair_good <- stair_paced %>%
  dplyr::filter(belt_quality == "good")
visual_good <- visual_ses1 %>%
  dplyr::filter(belt_quality == "good")

cat(sprintf("stair_paced (good belt): %d trials / %d participants\n",
            nrow(stair_good), dplyr::n_distinct(stair_good$id)))

sens1_A <- compare_sensitivity(
  model_A, stair_good,
  accuracy ~ direction_match + change + change_sq + (1 | id),
  family = binomial,
  label  = "H-physio-A (good belt)")

sens1_B <- compare_sensitivity(
  model_B, visual_good,
  accuracy ~ direction_match + change + change_sq + (1 | id),
  family = binomial,
  label  = "H-physio-B falsification (good belt)")

sens1_D <- compare_sensitivity(
  model_D, stair_good,
  arousal ~ direction_match + change + change_sq + (1 | id),
  label = "H-physio-D (good belt)")


# ── Sensitivity 2: Minimum slope trials (>= 25 of 32 non-NoChange) ───────────
cat("\n\n══════════════════════════════════════════════\n")
cat("SENSITIVITY 2 — Minimum slope trials (>= 25)\n")
cat("══════════════════════════════════════════════\n")

# Compute n_slope_trials per participant × session
slope_counts <- breathData %>%
  dplyr::filter(trial_type == "staircase", is_paced,
                change != 0, !is.na(ibi_slope)) %>%
  dplyr::group_by(id, ses) %>%
  dplyr::summarise(n_slope = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n_slope >= 25)

stair_min_slope <- stair_paced %>%
  dplyr::semi_join(slope_counts, by = c("id", "ses"))

# slope_counts was built on is_paced == TRUE, so Visual ses1 participants
# (who weren't pacing) are absent from it entirely — semi_join returns
# an empty tibble. Build a separate count for Visual ses1 without the
# is_paced filter, using the same n >= 25 threshold.
slope_counts_visual <- breathData %>%
  dplyr::filter(trial_type == "staircase",
                ses == 1,
                first_condition == "Visual",
                change != 0,
                !is.na(ibi_slope)) %>%
  dplyr::group_by(id, ses) %>%
  dplyr::summarise(n_slope = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n_slope >= 25)

visual_min_slope <- visual_ses1 %>%
  dplyr::semi_join(slope_counts_visual, by = c("id", "ses"))

cat(sprintf("stair_paced (n_slope >= 25): %d trials / %d participants\n",
            nrow(stair_min_slope), dplyr::n_distinct(stair_min_slope$id)))

sens2_A <- compare_sensitivity(
  model_A, stair_min_slope,
  accuracy ~ direction_match + change + change_sq + (1 | id),
  family = binomial,
  label  = "H-physio-A (>=25 slope trials)")

sens2_B <- compare_sensitivity(
  model_B, visual_min_slope,
  accuracy ~ direction_match + change + change_sq + (1 | id),
  family = binomial,
  label  = "H-physio-B falsification (>=25 slope trials)")


# ── Sensitivity 3: Belt quality and n_slope_trials as covariates ─────────────
cat("\n\n══════════════════════════════════════════════\n")
cat("SENSITIVITY 3 — Signal quality as covariate\n")
cat("══════════════════════════════════════════════\n")

# Join session-level slope count back to trial data
stair_with_qc <- stair_paced %>%
  dplyr::left_join(
    breathData %>%
      dplyr::filter(trial_type == "staircase", is_paced,
                    change != 0, !is.na(ibi_slope)) %>%
      dplyr::group_by(id, ses) %>%
      dplyr::summarise(n_slope_ses = dplyr::n(), .groups = "drop"),
    by = c("id", "ses")) %>%
  dplyr::mutate(
    belt_good      = as.integer(belt_quality == "good"),
    n_slope_ses_c  = scale(n_slope_ses, center = TRUE, scale = TRUE)[, 1]
  )

sens3_A <- glmer(
  accuracy ~ direction_match + change + change_sq +
             belt_good + n_slope_ses_c + (1 | id),
  data    = stair_with_qc,
  family  = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
cat("\nH-physio-A with belt quality and n_slope as covariates:\n")
print(summary(sens3_A)$coefficients)

sens3_D <- lmer(
  arousal ~ direction_match + change + change_sq +
            belt_good + n_slope_ses_c + (1 | id),
  data = stair_with_qc
)
cat("\nH-physio-D with belt quality and n_slope as covariates:\n")
print(summary(sens3_D)$coefficients)


# ── Sensitivity 4: Falsification specificity (z-test comparing β_Breath vs β_Visual) ──
cat("\n\n══════════════════════════════════════════════\n")
cat("SENSITIVITY 4 — Falsification specificity z-test\n")
cat("══════════════════════════════════════════════\n")

# Extract β and SE for direction_match from H-physio-A (paced) and B (Visual ses1)
coef_A <- summary(model_A)$coefficients["direction_match", ]
coef_B <- summary(model_B)$coefficients["direction_match", ]

beta_A <- coef_A["Estimate"]; se_A <- coef_A["Std. Error"]
beta_B <- coef_B["Estimate"]; se_B <- coef_B["Std. Error"]

z_diff <- (beta_A - beta_B) / sqrt(se_A^2 + se_B^2)
p_diff <- 2 * pnorm(-abs(z_diff))   # two-sided

cat(sprintf(
  "\nH-physio-A (paced) direction_match:   β = %.4f  SE = %.4f\n",
  beta_A, se_A))
cat(sprintf(
  "H-physio-B (Visual ses1) direction_match: β = %.4f  SE = %.4f\n",
  beta_B, se_B))
cat(sprintf(
  "Difference z-test: z = %.3f, p = %.4f\n", z_diff, p_diff))
cat("Significant difference supports specificity of adherence effect to paced sessions.\n")


# ── IBI drift covariate check ────────────────────────────────────────────────
cat("\n\n══════════════════════════════════════════════\n")
cat("IBI DRIFT COVARIATE — mean_ibi_slope as nuisance\n")
cat("══════════════════════════════════════════════\n")

# The positive mean IBI slope (0.041 s/breath, p = .0003) indicates a
# slight within-trial drift toward slower breathing independent of demand.
# Check whether controlling for mean_ibi_slope per session changes
# the direction_match effect on accuracy.

ibi_drift <- adhereSummary %>%
  dplyr::select(id, ses, mean_ibi_slope) %>%
  dplyr::rename(ibi_drift = mean_ibi_slope)

stair_drift <- stair_paced %>%
  dplyr::left_join(ibi_drift, by = c("id", "ses")) %>%
  dplyr::mutate(ibi_drift_c = scale(ibi_drift, center = TRUE, scale = TRUE)[, 1])

drift_A <- glmer(
  accuracy ~ direction_match + change + change_sq + ibi_drift_c + (1 | id),
  data    = stair_drift,
  family  = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
cat("\nH-physio-A controlling for session-level IBI drift:\n")
print(summary(drift_A)$coefficients)
cat(sprintf("direction_match β change: %.4f → %.4f\n",
            fixef(model_A)["direction_match"],
            fixef(drift_A)["direction_match"]))

# ══════════════════════════════════════════════════════════════════════════════
#  MODEL C ROBUSTNESS CHECKS
#
#  The ibi_slope_std × change interaction (β = 0.706) could partly reflect
#  that larger |Change| trials produce stronger IBI slopes simply because the
#  pacer moves more — making it easier to follow AND easier to detect.
#  Two checks:
#    RC1: Add |change| as a covariate alongside the interaction
#    RC2: Run the interaction separately by Salience condition
#         (Low = gradual change; High = abrupt change)
#         If the interaction holds in both salience conditions, it can't be
#         explained by overall signal strength alone.
# ══════════════════════════════════════════════════════════════════════════════

cat("\n\n══════════════════════════════════════════════\n")
cat("MODEL C ROBUSTNESS CHECKS\n")
cat("══════════════════════════════════════════════\n")

# ── RC1: Add |change| as additional covariate ─────────────────────────────────
cat("\n--- RC1: ibi_slope_std × change, controlling for |change| ---\n")
cat("Rationale: if larger |change| independently drives both better following\n")
cat("and better detection, the interaction should attenuate when |change| is\n")
cat("partialled out.\n\n")

stair_paced <- stair_paced %>%
  dplyr::mutate(abs_change = abs(change))

model_C_rc1 <- glmer(
  accuracy ~ ibi_slope_std_c * change + change_sq + abs_change + (1 | id),
  data    = stair_paced,
  family  = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
print(summary(model_C_rc1)$coefficients)

cat(sprintf(
  "\nInteraction β: original = %.4f  |  with |change| covariate = %.4f\n",
  fixef(model_C)["ibi_slope_std_c:change"],
  fixef(model_C_rc1)["ibi_slope_std_c:change"]))

cat("\nModel comparison (with vs without |change|):\n")
print(anova(model_C, model_C_rc1))


# ── RC2: Interaction by Salience condition ─────────────────────────────────────
cat("\n--- RC2: ibi_slope_std × change interaction by Salience ---\n")
cat("Rationale: Low salience = gradual change across 4 breaths;\n")
cat("High salience = abrupt change at breath 3.\n")
cat("If the interaction is mechanistic (not just signal strength),\n")
cat("it should hold within each salience condition.\n\n")

for (sal in c("Low", "High")) {
  sub <- stair_paced %>% dplyr::filter(salience == sal)
  cat(sprintf("Salience = %s (N trials = %d, N participants = %d):\n",
              sal, nrow(sub), dplyr::n_distinct(sub$id)))

  m_sal <- tryCatch(
    glmer(
      accuracy ~ ibi_slope_std_c * change + change_sq + (1 | id),
      data    = sub,
      family  = binomial,
      control = glmerControl(optimizer = "bobyqa")
    ),
    error = function(e) {
      message(sprintf("  Model failed: %s", e$message)); NULL
    }
  )

  if (!is.null(m_sal)) {
    coefs <- summary(m_sal)$coefficients
    cat(sprintf(
      "  ibi_slope_std_c:change  β = %.4f  SE = %.4f  z = %.3f  p = %.4f\n\n",
      coefs["ibi_slope_std_c:change", "Estimate"],
      coefs["ibi_slope_std_c:change", "Std. Error"],
      coefs["ibi_slope_std_c:change", "z value"],
      coefs["ibi_slope_std_c:change", "Pr(>|z|)"]))
  }
}


# ── RC3: Three-way interaction ibi_slope_std × change × Salience ─────────────
cat("\n--- RC3: Three-way ibi_slope_std × change × Salience ---\n")
cat("Directly tests whether the interaction differs between salience conditions.\n\n")

stair_paced <- stair_paced %>%
  dplyr::mutate(salience_c = ifelse(salience == "High", 0.5, -0.5))  # sum coding

model_C_3way <- glmer(
  accuracy ~ ibi_slope_std_c * change * salience_c + change_sq + (1 | id),
  data    = stair_paced,
  family  = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
print(summary(model_C_3way)$coefficients)

cat("\nModel comparison (2-way vs 3-way):\n")
print(anova(model_C, model_C_3way))
