# ============================================================
#  Intero2025_TrimContaminatedPhysio.R
#
#  Handles physio files that contain a previous participant's
#  recording prepended to the actual session.
#
#  ROOT CAUSE: The BioPac template was not cleared between
#  sessions on certain dates, causing the previous recording
#  to be included at the start of the next file.
#
#  AFFECTED PARTICIPANTS (confirmed or suspected):
#
#  Trigger-based trim (abcode = 3: one extra session trigger):
#    14359  R  inherits from 15571 (last R before Oct 7)
#    15313  R  inherits from 14359
#    15961  R  inherits from 15313
#    15367  R  inherits from unknown (1 trigger found earlier;
#               re-diagnosed here for consistency)
#
#  Signal-match trim (abcode = 2, previous had no trigger):
#    15091  R  inherits from 14338  (RA note confirms)
#    14581  L  inherits from 16456  (last L-side before Oct 7)
#    14101  L  inherits from 14581
#
#  Skipped / no action:
#    16510  no physio file found
#    12982  RA note: ".acq is fine"
#    Large paired files (9238, 15724, 14764, 13771, 12598,
#      13531): legitimate simultaneous L+R recordings, not
#      contaminated. checkParticipantFiles already handles
#      channel separation correctly.
#
#  Needs inspection (abcode = 3 but no contamination note):
#    15025, 9892  — diagnose only, no auto-trim
#    10348/7972   — 4 triggers, 116-min recording, inspect
#
#  HOW TO USE:
#    1. Source this file (paths set below).
#    2. Run STEP 1 (diagnostics) — always read-only.
#    3. Review dry-run output and diagnostic plots in
#       contamination_qc/.
#    4. Change dry_run = FALSE in STEP 3 for each participant
#       only after you are satisfied with the dry run.
#
#  Outputs:
#    contamination_qc/   — diagnostic plots (signal-match)
#    rds/<id>.rds        — updated: trimmed signal, adjusted
#                          start_indices, $contamination_note
# ============================================================

# Set Up ---------
## Load libraries ---------
packages <- c("ggplot2", "patchwork", "signal", "reticulate")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages)
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)


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

source(file.path(analysisPath, "breath_pipeline.R"))

reticulate::py_install("bioread")
bioread <- reticulate::import("bioread")

rds_dir  <- file.path(resultsPath, "rds")
plot_dir <- file.path(resultsPath, "contamination_qc")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# File lookup (from belt_quality_screen.csv)
physio_files <- list(
  "14359" = "L14581.R14359.acq",
  "15313" = "L14101.R15313.acq",
  "15961" = "l16510.R15961.acq",   # note lowercase l in filename
  "15367" = "L00000.R15367.acq",
  "14581" = "L14581.R14359.acq",
  "14101" = "L14101.R15313.acq",
  "15091" = "L0000.R15091.acq",
  "15025" = "L15199.R15025.acq",
  "9892"  = "L9892.R3472.acq",
  "10348" = "L7972.R10348.acq",
  "7972"  = "L7972.R10348.acq"
)


# ============================================================
#  HELPER: diagnose_acq_triggers
#  Read the raw .acq trigger channel and report all events.
#  Always run before trimming any file.
# ============================================================
diagnose_acq_triggers <- function(acq_path, bioread,
                                  trigger_ch_idx   = 13L,
                                  participant_label = "") {
  message(sprintf("\n=== Trigger diagnostics: %s ===", participant_label))

  dat <- tryCatch(bioread$read_file(acq_path),
                  error = function(e) { message("  [ERROR] ", e$message); NULL })
  if (is.null(dat)) return(invisible(NULL))

  sr     <- as.numeric(dat$samples_per_second)
  n_ch   <- length(dat$channels)
  n_samp <- length(dat$channels[[1L]]$data)
  dur_s  <- n_samp / sr
  message(sprintf("  SR: %.0fHz | Channels: %d | Duration: %.1fs (%.1f min)",
                  sr, n_ch, dur_s, dur_s / 60))

  if (trigger_ch_idx > n_ch) {
    message(sprintf("  [WARN] Trigger channel %d not found", trigger_ch_idx))
    return(invisible(NULL))
  }

  triggers <- as.numeric(dat$channels[[trigger_ch_idx]]$data)
  changes  <- which(c(TRUE, diff(triggers) != 0) & triggers != 0 & triggers < 200)
  trig_df  <- data.frame(
    event    = seq_along(changes),
    sample   = changes,
    time_s   = round(changes / sr, 3),
    time_min = round(changes / sr / 60, 2),
    value    = triggers[changes]
  )
  message(sprintf("  Total trigger events: %d (all values)", nrow(trig_df)))

  ses_starts <- trig_df[trig_df$value %in% c(1L, 16L), ]
  message(sprintf("  Session-start triggers (value 1 or 16): %d", nrow(ses_starts)))
  if (nrow(ses_starts) > 0) print(ses_starts) else message("  None found.")

  invisible(trig_df)
}


# ============================================================
#  HELPER: trim_rds_by_trigger
#
#  For participants whose file has one extra session-start
#  trigger due to template contamination (abcode = 3).
#
#  Trigger assignment by count of session-start triggers
#  of the expected value (1 for R-side, 16 for L-side):
#
#    1 trigger (Scenario A):
#      The trigger belongs to THIS participant's ses1.
#      Trim to it; no ses2 trigger exists.
#      new start_indices = [0]
#
#    2 triggers (Scenario B):
#      Trigger 1 = previous participant's ses1.
#      Trigger 2 = this participant's ses1.
#      No ses2 trigger.
#      new start_indices = [ses1_ds_from_new_t0]
#
#    3 triggers (Scenario C — Oct 7 cascade pattern):
#      Trigger 1 = previous participant's ses1.
#      Trigger 2 = this participant's ses1.
#      Trigger 3 = this participant's ses2.
#      new start_indices = [ses1_ds, ses2_ds] from new t=0.
#
#    4+ triggers: printed for manual inspection, no action.
#
#  In all cases the recording is trimmed to start at the
#  1st trigger so the contaminated pre-recording is removed.
# ============================================================
trim_rds_by_trigger <- function(participant_id,
                                acq_path,
                                bioread,
                                rds_dir,
                                side           = "R",
                                trigger_val_R  = 1L,
                                trigger_val_L  = 16L,
                                trigger_ch_idx = 13L,
                                dry_run        = TRUE) {

  message(sprintf("\n=== Trigger-based trim: %s | side=%s | dry_run=%s ===",
                  participant_id, side, dry_run))

  dat       <- bioread$read_file(acq_path)
  native_sr <- as.numeric(dat$samples_per_second)
  triggers  <- as.numeric(dat$channels[[trigger_ch_idx]]$data)
  changes   <- which(c(TRUE, diff(triggers) != 0) & triggers != 0 & triggers < 200)
  trig_vals <- triggers[changes]

  expected_val <- if (side == "R") trigger_val_R else trigger_val_L
  ses_idx      <- which(trig_vals == expected_val)
  n_ses        <- length(ses_idx)

  rds_obj <- readRDS(file.path(rds_dir, paste0(participant_id, ".rds")))
  ds_fs   <- rds_obj$fs
  ds_q    <- floor(native_sr / ds_fs)

  message(sprintf("  Native SR: %.0fHz | ds_fs: %.0fHz | Session-start triggers (val=%d): %d",
                  native_sr, ds_fs, expected_val, n_ses))

  if (n_ses == 0) {
    message("  No session-start triggers found.")
    return(invisible(NULL))
  }

  # Always trim to the 1st trigger (start of contaminated block)
  trim_native <- changes[ses_idx[1]]
  trim_ds     <- round(trim_native / ds_q)
  trim_s      <- trim_native / native_sr

  if (n_ses == 1) {
    # The single trigger at t≈173.2s (sample 346448) is the R-side ghost
    # baked into the template — it belongs to a previous participant.
    # After trimming, no real session triggers remain.  Set start_indices
    # to NULL so alignment recovery finds session starts from the breathing
    # pattern rather than using a wrong onset.
    new_starts <- NULL
    scenario   <- "A (1 ghost trigger only — start_indices NULL, alignment recovery required)"

  } else if (n_ses == 2) {
    ses1_ds    <- round(changes[ses_idx[2]] / ds_q) - trim_ds
    new_starts <- ses1_ds
    scenario   <- sprintf("B (2 triggers: prev ses1 at trim point, this ses1 at new t+%.1fs)",
                          ses1_ds / ds_fs)

  } else if (n_ses == 3) {
    ses1_ds    <- round(changes[ses_idx[2]] / ds_q) - trim_ds
    ses2_ds    <- round(changes[ses_idx[3]] / ds_q) - trim_ds
    new_starts <- c(ses1_ds, ses2_ds)
    scenario   <- sprintf(
      "C (3 triggers: prev at trim, ses1 at +%.1fs, ses2 at +%.1fs)",
      ses1_ds / ds_fs, ses2_ds / ds_fs)

  } else if (n_ses == 4) {
    # Scenario D: two full previous-participant sessions prepended.
    # Trigger 1 = previous ses1, trigger 2 = previous ses2,
    # trigger 3 = this participant's ses1, trigger 4 = this participant's ses2.
    # Trim to trigger 1; real ses1/ses2 are triggers 3 and 4.
    ses1_ds    <- round(changes[ses_idx[3]] / ds_q) - trim_ds
    ses2_ds    <- round(changes[ses_idx[4]] / ds_q) - trim_ds
    new_starts <- c(ses1_ds, ses2_ds)
    scenario   <- sprintf(
      "D (4 triggers: 2 prev + ses1 at +%.1fs + ses2 at +%.1fs)",
      ses1_ds / ds_fs, ses2_ds / ds_fs)

  } else {
    message(sprintf("  %d session-start triggers — no handler, inspect manually.", n_ses))
    print(data.frame(sample  = changes[ses_idx],
                     time_s  = round(changes[ses_idx] / native_sr, 2),
                     value   = trig_vals[ses_idx]))
    return(invisible(NULL))
  }

  remain_s <- (length(rds_obj$signal) - trim_ds) / ds_fs
  message(sprintf("  Scenario: %s", scenario))
  message(sprintf("  Trim at native sample %d (t=%.2fs)", trim_native, trim_s))
  message(sprintf("  Remaining signal: %.1fs (%.1f min)", remain_s, remain_s / 60))
  message(sprintf("  New start_indices (ds samples): %s",
                  paste(new_starts, collapse = ", ")))

  if (dry_run) {
    message("  [DRY RUN] No files modified.")
    return(invisible(list(scenario = scenario, trim_s = trim_s,
                          trim_ds = trim_ds, new_starts = new_starts,
                          remain_s = remain_s)))
  }

  rds_obj$signal             <- rds_obj$signal[(trim_ds + 1L):length(rds_obj$signal)]
  rds_obj$trim_offset_s      <- trim_s
  rds_obj$trim_offset_samp   <- trim_ds
  rds_obj$start_indices      <- new_starts
  rds_obj$contamination_note <- sprintf(
    "Trigger trim: %.1fs removed. %s.", trim_s, scenario)

  saveRDS(rds_obj, file.path(rds_dir, paste0(participant_id, ".rds")))
  message(sprintf("  [SAVED] %s.rds", participant_id))
  invisible(rds_obj)
}


# ============================================================
#  HELPERS: signal-match trim
#  For abcode=2 participants where the previous session left
#  no trigger — find transition via sliding cross-correlation.
# ============================================================
find_signal_match_offset <- function(long_signal, short_signal, ds_fs,
                                     search_window_s = 300,
                                     validation_r    = 0.90,
                                     segment_s       = 60) {
  seg_samp  <- min(round(segment_s * ds_fs), length(short_signal))
  short_seg <- short_signal[seq_len(seg_samp)]
  short_seg <- (short_seg - mean(short_seg)) / sd(short_seg)

  search_samp <- min(round(search_window_s * ds_fs),
                     length(long_signal) - seg_samp)
  message(sprintf("  Sliding correlation: %d samples (%.0f s) search window...",
                  search_samp, search_samp / ds_fs))

  cors <- vapply(seq_len(search_samp), function(i) {
    w <- long_signal[i:(i + seg_samp - 1L)]
    s <- sd(w)
    if (is.na(s) || s < 1e-10) return(NA_real_)
    cor((w - mean(w)) / s, short_seg, use = "complete.obs")
  }, numeric(1))

  best_idx <- which.max(cors)
  best_r   <- cors[best_idx]
  best_t   <- best_idx / ds_fs
  message(sprintf("  Best match: sample %d (t=%.2fs), r=%.4f", best_idx, best_t, best_r))

  if (best_r < validation_r)
    warning(sprintf("[WARN] r=%.3f below threshold %.2f — inspect plot carefully.",
                    best_r, validation_r))
  else
    message(sprintf("  Match validated (r=%.4f)", best_r))

  list(start_sample = best_idx, start_s = best_t, best_r = best_r,
       cors = cors, ds_fs = ds_fs)
}

plot_match_diagnostic <- function(match_result, long_signal, short_signal,
                                  ds_fs, participant_label, save_path = NULL) {
  mm <- function(x) {
    r <- range(x, na.rm = TRUE)
    if (diff(r) < 1e-10) return(rep(0, length(x)))
    2 * (x - r[1]) / diff(r) - 1
  }
  n_show  <- min(round(60 * ds_fs), length(short_signal),
                 length(long_signal) - match_result$start_sample + 1L)
  t_show  <- seq_len(n_show) / ds_fs

  df_cor <- data.frame(t_s = seq_along(match_result$cors) / ds_fs,
                       r   = match_result$cors)
  p1 <- ggplot(df_cor, aes(x = t_s, y = r)) +
    geom_line(colour = "#2c7be5", linewidth = 0.5) +
    geom_vline(xintercept = match_result$start_s,
               colour = "darkorange", linetype = "dashed") +
    annotate("text", x = match_result$start_s,
             y = max(df_cor$r, na.rm = TRUE),
             label = sprintf("t=%.2fs\nr=%.4f",
                             match_result$start_s, match_result$best_r),
             hjust = -0.1, vjust = 1, size = 3, colour = "darkorange") +
    labs(title = sprintf("%s — sliding correlation landscape", participant_label),
         x = "Candidate offset (s)", y = "Pearson r") +
    theme_minimal(base_size = 10)

  obs_seg <- long_signal[match_result$start_sample:(match_result$start_sample + n_show - 1L)]
  ref_seg <- short_signal[seq_len(n_show)]
  df_ov   <- data.frame(
    t_s    = rep(t_show, 2),
    signal = c(mm(obs_seg), mm(ref_seg)),
    source = rep(c(sprintf("%s (matched)", participant_label),
                   "Reference (previous participant)"), each = n_show)
  )
  p2 <- ggplot(df_ov, aes(x = t_s, y = signal, colour = source)) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = c("#2c7be5", "#e63946")) +
    labs(title = "Signal overlay at match (first 60 s)",
         x = "Time (s)", y = "Normalised signal", colour = NULL) +
    theme_minimal(base_size = 10) + theme(legend.position = "top")

  combined <- patchwork::wrap_plots(p1, p2, ncol = 1)
  if (!is.null(save_path)) {
    ggsave(save_path, combined, width = 12, height = 8, dpi = 120)
    message(sprintf("  Plot saved: %s", basename(save_path)))
  } else print(combined)
  invisible(combined)
}

trim_rds_by_signal_match <- function(contaminated_id, previous_id,
                                     rds_dir, plot_dir,
                                     search_window_s = 300,
                                     validation_r    = 0.90,
                                     dry_run         = TRUE) {
  message(sprintf("\n=== Signal-match trim: %s (prev: %s) | dry_run=%s ===",
                  contaminated_id, previous_id, dry_run))

  rds_cont <- readRDS(file.path(rds_dir, paste0(contaminated_id, ".rds")))
  rds_prev <- readRDS(file.path(rds_dir, paste0(previous_id,     ".rds")))
  ds_fs    <- rds_cont$fs
  if (rds_prev$fs != ds_fs) stop("Sampling rates differ between RDS files.")

  message(sprintf("  Contaminated: %.1fs (%.1f min) | Previous: %.1fs (%.1f min)",
                  length(rds_cont$signal) / ds_fs, length(rds_cont$signal) / ds_fs / 60,
                  length(rds_prev$signal) / ds_fs, length(rds_prev$signal) / ds_fs / 60))

  match_r <- find_signal_match_offset(
    rds_cont$signal, rds_prev$signal, ds_fs,
    search_window_s = search_window_s, validation_r = validation_r)

  # Trim point = where previous participant's signal ends
  trim_samp <- match_r$start_sample + length(rds_prev$signal)
  trim_s    <- trim_samp / ds_fs
  remain_s  <- (length(rds_cont$signal) - trim_samp) / ds_fs
  message(sprintf("  Previous signal ends at sample %d (t=%.2fs / %.1f min)",
                  trim_samp, trim_s, trim_s / 60))
  message(sprintf("  Remaining after trim: %.1fs (%.1f min)", remain_s, remain_s / 60))

  # Adjust any pre-existing start_indices relative to new t=0
  new_starts <- if (!is.null(rds_cont$start_indices))
    rds_cont$start_indices[rds_cont$start_indices > trim_samp] - trim_samp
  else NULL
  message(sprintf("  Adjusted start_indices: %s",
                  if (length(new_starts)) paste(new_starts, collapse=", ") else "none"))

  plot_match_diagnostic(
    match_r, rds_cont$signal, rds_prev$signal, ds_fs,
    participant_label = as.character(contaminated_id),
    save_path = file.path(plot_dir,
                          sprintf("P%s_contamination_match.png", contaminated_id))
  )

  if (dry_run) {
    message("  [DRY RUN] No files modified.")
    return(invisible(list(match = match_r, trim_samp = trim_samp,
                          trim_s = trim_s, remain_s = remain_s,
                          new_starts = new_starts)))
  }

  rds_cont$signal             <- rds_cont$signal[(trim_samp + 1L):length(rds_cont$signal)]
  rds_cont$trim_offset_s      <- trim_s
  rds_cont$trim_offset_samp   <- trim_samp
  rds_cont$start_indices      <- new_starts
  rds_cont$contamination_note <- sprintf(
    "Signal-match trim: %.1fs removed (prev=%s, r=%.4f)",
    trim_s, previous_id, match_r$best_r)

  saveRDS(rds_cont, file.path(rds_dir, paste0(contaminated_id, ".rds")))
  message(sprintf("  [SAVED] %s.rds", contaminated_id))
  invisible(rds_cont)
}


# ============================================================
#  STEP 1 — DIAGNOSTICS (read-only, always run first)
# ============================================================
message("\n\n==================== STEP 1: DIAGNOSTICS ====================\n")

# Oct 7 R-side cascade: expect 3 session-start triggers each
for (pid in c("14359", "15313", "15961"))
  diagnose_acq_triggers(file.path(physioPath, physio_files[[pid]]),
                        bioread, participant_label = pid)

# 15367: check whether raw file has 1 or 2 session-start triggers
diagnose_acq_triggers(file.path(physioPath, physio_files[["15367"]]),
                      bioread, participant_label = "15367")

# Inspect-only (no contamination note but unusual trigger count)
for (pid in c("15025", "9892", "10348", "7972"))
  diagnose_acq_triggers(file.path(physioPath, physio_files[[pid]]),
                        bioread, participant_label = pid)

# Interpretation of inspect-only cases from trigger table:
#
#  15025, 9892: Normal duration (~59 min), no RA contamination note.
#    Extra trigger is likely cable noise (~80s after real ses1).
#    No trimming needed — their abcode=3 is a false alarm.
#
#  10348 (R-side) / 7972 (L-side): Contaminated (116.8 min, 4 triggers).
#    R-side (val=1) triggers: 200.1s*, 2132.8s*, 3696.4s, 5800.7s
#    L-side (val=16) triggers: 122.5s*, 326.4s*, 2411.3s, 5632.6s
#    * = previous participant's triggers; real sessions start at trigger 3/4.
#    Trim 10348 to t=200.1s (1st val=1); real ses1=3696.4s, ses2=5800.7s.
#    Trim 7972  to t=122.5s (1st val=16); real ses1=2411.3s, ses2=5632.6s.
#    Handled explicitly in Step 2/3 below.


# ============================================================
#  STEP 2 — DRY RUNS (review output before proceeding)
# ============================================================
message("\n\n==================== STEP 2: DRY RUNS ====================\n")

# --- Trigger-based (Oct 7 R-side + 15367) ---
for (pid in c("14359", "15313", "15961", "15367")) {
  trim_rds_by_trigger(
    participant_id = pid,
    acq_path       = file.path(physioPath, physio_files[[pid]]),
    bioread        = bioread,
    rds_dir        = rds_dir,
    side           = "R",
    dry_run        = TRUE
  )
}

# --- 10348 (R-side) and 7972 (L-side) — confirmed contamination ---
# Both have 4 session-start triggers; the function falls into the manual
# inspection branch.  Handled here by explicit side assignment so the
# function picks the correct trigger value for each.
trim_rds_by_trigger(
  participant_id = "10348",
  acq_path       = file.path(physioPath, physio_files[["10348"]]),
  bioread        = bioread, rds_dir = rds_dir,
  side    = "R",   # val=1 triggers: ghost at 200.1s, ghost ses2 2132.8s, real ses1 3696.4s, real ses2 5800.7s
  dry_run = TRUE
)
trim_rds_by_trigger(
  participant_id = "7972",
  acq_path       = file.path(physioPath, physio_files[["7972"]]),
  bioread        = bioread, rds_dir = rds_dir,
  side    = "L",   # val=16 triggers: ghost ses1 122.5s, ghost ses2 326.4s, real ses1 2411.3s, real ses2 5632.6s
  dry_run = TRUE
)

# --- Signal-match (L-side, no trigger; trim IN ORDER) ---
# 15091: previous = 14338 (RA confirmed, r=1.000 at sample 1)
trim_rds_by_signal_match("15091", "14338", rds_dir, plot_dir,
                         search_window_s = 300, dry_run = TRUE)

# 14581 and 14101: both L-side (val=16) with abcode=2 — two triggers each.
# 16456 (last L-session before Oct 7) has no physio file so signal-match
# is unavailable.  Trigger-based trim applies instead (Scenario B):
# trigger 1 = previous participant's ses1, trigger 2 = this participant's ses1.
trim_rds_by_trigger(
  participant_id = "14581",
  acq_path       = file.path(physioPath, physio_files[["14581"]]),
  bioread        = bioread, rds_dir = rds_dir,
  side    = "L",
  dry_run = TRUE
)
trim_rds_by_trigger(
  participant_id = "14101",
  acq_path       = file.path(physioPath, physio_files[["14101"]]),
  bioread        = bioread, rds_dir = rds_dir,
  side    = "L",
  dry_run = TRUE
)


# ============================================================
#  STEP 3 — LIVE TRIMS
#  Uncomment each block only after reviewing dry-run output.
#  For signal-match, trim IN ORDER (14581 before 14101).
# ============================================================
message("\n\n==================== STEP 3: LIVE TRIMS ====================\n")
message("All live trims are commented out. Uncomment after reviewing dry runs.\n")

# Trigger-based
for (pid in c("14359", "15313", "15961", "15367")) {
  trim_rds_by_trigger(
    participant_id = pid,
    acq_path       = file.path(physioPath, physio_files[[pid]]),
    bioread        = bioread, rds_dir = rds_dir,
    side = "R", dry_run = FALSE
  )
}
trim_rds_by_trigger("14581",
                    acq_path=file.path(physioPath,physio_files[["14581"]]),
                    bioread=bioread, rds_dir=rds_dir, side="L", dry_run=FALSE)
trim_rds_by_trigger("14101",
                    acq_path=file.path(physioPath,physio_files[["14101"]]),
                    bioread=bioread, rds_dir=rds_dir, side="L", dry_run=FALSE)
trim_rds_by_trigger("10348",
                    acq_path=file.path(physioPath,physio_files[["10348"]]),
                    bioread=bioread, rds_dir=rds_dir, side="R", dry_run=FALSE)
trim_rds_by_trigger("7972",
                    acq_path=file.path(physioPath,physio_files[["7972"]]),
                    bioread=bioread, rds_dir=rds_dir, side="L", dry_run=FALSE)

# Signal-match
trim_rds_by_signal_match("15091", "14338", rds_dir, plot_dir,
                         search_window_s=300, dry_run=FALSE)