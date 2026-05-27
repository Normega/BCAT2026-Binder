# =============================================================
#  Intero2025_AlignmentRecovery.R
#
#  Recovers the alignment between task log data and physio
#  recording for participants where trigger signals are absent
#  or unreliable.
#
#  Strategy:
#    The pipeline has already detected every breath trough in
#    the recording.  For any proposed run-start offset t, we
#    can synthesise exactly where every pacer trough *should*
#    have occurred (from the log data).  We score the alignment
#    by the mean absolute error between expected and nearest
#    observed troughs.  The offset that minimises this MAE is
#    the recovered alignment.
#
#  Two-stage search per run:
#    1. Coarse pass (1 s steps) over the full valid range
#    2. Fine pass  (0.1 s steps) within ±fine_window_s of best
#
#  Two-run handling:
#    firstCondition == "Breath" → ses1 AND ses2 are paced:
#      sequential search with constraint t2 > t1 + run1_span + min_gap
#      followed by joint 2-D refinement
#    firstCondition == "Visual" → only ses2 is paced:
#      t2 lower bound set from trigger_onset_ses1_s when available,
#      otherwise from visual_session_span_s + min_gap_s
#      (ses1 occupies the recording but produces no expected troughs)
#
#  Changes from original version:
#    - score_trough_alignment vectorised with outer() + apply(,2,min)
#      instead of vapply loop — ~10x faster for typical trough counts
#    - compute_trial_durations exponent uses (NUMBREATHS - 1) instead
#      of hardcoded 3 — safe if NUMBREATHS ever changes from 4
#    - Visual ses1 lower bound sharpened: uses trigger_onset_ses1_s
#      when available so the search range is tighter and faster
#    - confidence renamed to ambiguity throughout (higher = worse;
#      ambiguity > 0.7 flags a flat / unresolved cost landscape)
#    - belt_quality parameter: returns early with a note if "unusable"
#      so the cost of a full search is never paid on dead signals
#
#  Dependencies:
#    breath_pipeline.R   — for run_pipeline() / trough_times
#    Intero2025_RespirationFunctions.R — for getExpectedSignal()
# =============================================================


# -------------------------------------------------------------
#  compute_trial_durations
#
#  Replicates the breath-duration logic from getExpectedSignal()
#  without needing the signal vector.  Returns a numeric vector
#  of length NUMBREATHS (one duration per breath in the trial).
#
#  Arguments:
#    Salience   — "Low" or "High"
#    Change     — numeric change factor (0 = no change)
#    STARTDUR   — base breath duration in seconds (default 4)
#    NUMBREATHS — breaths per trial (default 4); used to compute
#                 the geometric step for Low salience trials.
#                 Must match the value used in getExpectedSignal().
# -------------------------------------------------------------
compute_trial_durations <- function(Salience, Change,
                                    STARTDUR   = 4,
                                    NUMBREATHS = 4) {

  if (any(is.na(Salience)) || any(is.na(Change)) ||
      length(Salience) == 0 || length(Change) == 0)
    return(rep(NA_real_, NUMBREATHS))

  Salience <- Salience[1]
  Change   <- as.numeric(Change[1])

  changeVal    <- 1 + Change
  # FIX: was hardcoded as ^ (1/3); now uses (NUMBREATHS - 1) so the
  # function stays correct if NUMBREATHS is ever changed from 4.
  changeAmount <- changeVal ^ (1 / (NUMBREATHS - 1))

  if (Change == 0)
    return(rep(STARTDUR, NUMBREATHS))

  if (Salience == "Low") {
    vapply(seq_len(NUMBREATHS), function(k)
      STARTDUR * changeAmount ^ (k - 1L), numeric(1))
  } else {
    # High salience: change applied abruptly between breath 2 and 3
    c(rep(STARTDUR, 2L),
      rep(STARTDUR * (1 + Change), NUMBREATHS - 2L))
  }
}


# -------------------------------------------------------------
#  normalise_test_data
#
#  Converts test-phase trial data (from testFile.csv) to the
#  same column format as staircase data so the two can be
#  combined and passed to synthesise_expected_troughs.
#
#  Test file uses:
#    TestTrial.started / TestTrial.stopped  (timing)
#    Salience: integer 0 (Low) / 1 (High)
#    Direction: integer -1 (Faster) / 0 (NoChange) / 1 (Slower)
#    level: positive numeric (magnitude of change)
#
#  Output columns match staircase format:
#    trial.started, trial.stopped, Salience ("Low"/"High"), Change
# -------------------------------------------------------------
normalise_test_data <- function(test_data) {
  if (is.null(test_data) || nrow(test_data) == 0) return(NULL)

  data.frame(
    trial.started = test_data$TestTrial.started,
    trial.stopped = test_data$TestTrial.stopped,
    Salience      = ifelse(test_data$Salience == 1, "High", "Low"),
    Change        = test_data$Direction * test_data$level,
    stringsAsFactors = FALSE
  )
}


# -------------------------------------------------------------
#  synthesise_expected_troughs
#
#  Given a data frame of trials and a proposed recording-
#  absolute offset for the start of that run, returns the
#  recording-absolute times of all expected pacer troughs.
#
#  Each trial contributes (NUMBREATHS + 1) trough times:
#    offset + trial.started + c(0, d1, d1+d2, ..., sum(d))
#
#  Arguments:
#    trial_data  — data.frame with columns:
#                    trial.started (s, relative to run start)
#                    Salience ("Low"/"High"), Change (signed numeric)
#    test_data   — optional data.frame of test-phase trials in the
#                    SAME normalised format (i.e. already passed
#                    through normalise_test_data).  When supplied,
#                    test trough times are appended.  NULL = ignored.
#    offset_s    — proposed recording-absolute run start (s)
#    STARTDUR    — base breath duration (default 4)
#    NUMBREATHS  — breaths per trial (default 4)
#
#  Returns a sorted numeric vector of expected trough times
#  in recording-absolute seconds.
# -------------------------------------------------------------
synthesise_expected_troughs <- function(trial_data, offset_s,
                                        STARTDUR   = 4,
                                        NUMBREATHS = 4,
                                        test_data  = NULL) {

  # Combine staircase and test-phase rows (test_data already normalised)
  combined <- if (!is.null(test_data) && nrow(test_data) > 0)
    rbind(trial_data[, c("trial.started", "Salience", "Change")],
          test_data[ , c("trial.started", "Salience", "Change")])
  else
    trial_data

  all_troughs <- unlist(lapply(seq_len(nrow(combined)), function(i) {
    row  <- combined[i, ]
    durs <- compute_trial_durations(row$Salience, row$Change,
                                    STARTDUR, NUMBREATHS)
    if (any(is.na(durs))) return(NULL)
    offset_s + row$trial.started + c(0, cumsum(durs))
  }))

  sort(as.numeric(all_troughs))
}


# -------------------------------------------------------------
#  score_trough_alignment
#
#  For each expected trough, find the nearest observed trough
#  and compute absolute timing error.  Expected troughs with
#  no observed trough within max_match_dist_s are penalised
#  at max_match_dist_s (smooth penalty rather than Inf, keeps
#  the cost landscape well-behaved).
#
#  Returns mean cost across ALL expected troughs (so the score
#  is comparable across different run lengths).
#
#  Implementation note:
#    Uses outer() to build the full pairwise distance matrix
#    then apply(, 2, min) for column minima — vectorised and
#    substantially faster than the previous vapply loop over
#    expected troughs for the large number of offset evaluations
#    made during the coarse search pass.
#
#  Arguments:
#    obs_trough_times — observed trough times (recording-abs s)
#    exp_trough_times — expected trough times from synthesise_*
#    max_match_dist_s — penalty cap and match threshold (default 2)
# -------------------------------------------------------------
score_trough_alignment <- function(obs_trough_times,
                                   exp_trough_times,
                                   max_match_dist_s = 2) {

  if (length(exp_trough_times) == 0) return(NA_real_)
  if (length(obs_trough_times) == 0) return(max_match_dist_s)

  # Full pairwise absolute-distance matrix: rows = observed, cols = expected.
  # apply(, 2, min) gives the nearest-observed distance for each expected trough.
  dist_mat <- abs(outer(obs_trough_times, exp_trough_times, "-"))
  nearest  <- apply(dist_mat, 2L, min)
  mean(pmin(nearest, max_match_dist_s), na.rm = TRUE)
}


# -------------------------------------------------------------
#  .two_pass_search
#
#  Internal helper.  Runs coarse + fine search for a single
#  run's offset given a fixed set of observed troughs and
#  a function that synthesises expected troughs at a given
#  offset.
#
#  Arguments:
#    obs_troughs      — observed trough times (recording-abs s)
#    synth_fn         — function(offset_s) → expected trough times
#    search_min_s     — lower bound of valid offset range
#    search_max_s     — upper bound of valid offset range
#    coarse_step_s    — step size for coarse pass (default 1)
#    fine_step_s      — step size for fine pass   (default 0.1)
#    fine_window_s    — half-width of fine search  (default 5)
#    max_match_dist_s — passed to score_trough_alignment
#
#  Returns list:
#    best_offset      — best offset in seconds
#    best_mae         — MAE at best offset
#    coarse_landscape — data.frame(offset, mae) for coarse pass
#    fine_landscape   — data.frame(offset, mae) for fine pass
#    ambiguity        — best_mae / median(coarse_mae):
#                       higher → flatter landscape → less confident
#                       Values near 0: clear unique minimum.
#                       Values near 1: flat landscape, ambiguous.
# -------------------------------------------------------------
.two_pass_search <- function(obs_troughs, synth_fn,
                             search_min_s, search_max_s,
                             coarse_step_s    = 1,
                             fine_step_s      = 0.1,
                             fine_window_s    = 5,
                             max_match_dist_s = 2) {

  # ---- Coarse pass ----
  coarse_offsets <- seq(search_min_s, search_max_s, by = coarse_step_s)

  if (length(coarse_offsets) == 0) {
    message("  [WARNING] Search range is empty — check recording duration and run spans.")
    return(NULL)
  }

  n_coarse    <- length(coarse_offsets)
  coarse_maes <- numeric(n_coarse)
  print_every <- max(1L, floor(n_coarse / 20))

  for (i in seq_len(n_coarse)) {
    coarse_maes[i] <- score_trough_alignment(
      obs_troughs, synth_fn(coarse_offsets[i]), max_match_dist_s)
    if (i %% print_every == 0 || i == n_coarse)
      cat(sprintf("\r    Coarse: %d/%d  (%.0f%%)  best MAE so far: %.3f s   ",
                  i, n_coarse, 100 * i / n_coarse,
                  min(coarse_maes[1:i], na.rm = TRUE)))
  }
  cat("\n")

  coarse_landscape <- data.frame(offset = coarse_offsets, mae = coarse_maes)
  best_coarse      <- coarse_offsets[which.min(coarse_maes)]

  # ---- Fine pass ----
  fine_min     <- max(search_min_s, best_coarse - fine_window_s)
  fine_max     <- min(search_max_s, best_coarse + fine_window_s)
  fine_offsets <- seq(fine_min, fine_max, by = fine_step_s)

  fine_maes <- vapply(fine_offsets, function(t)
    score_trough_alignment(obs_troughs, synth_fn(t), max_match_dist_s),
    numeric(1))

  fine_landscape <- data.frame(offset = fine_offsets, mae = fine_maes)
  best_fine      <- fine_offsets[which.min(fine_maes)]
  best_mae       <- min(fine_maes)

  # Ambiguity: ratio of best MAE to median of coarse landscape.
  # Higher value = flatter landscape = less certain alignment.
  # Values near 0 = clear unique minimum (low ambiguity).
  # Values near 1 = flat landscape     (high ambiguity).
  med_mae   <- median(coarse_maes, na.rm = TRUE)
  ambiguity <- if (!is.na(med_mae) && med_mae > 0)
    best_mae / med_mae else NA_real_

  list(
    best_offset      = best_fine,
    best_mae         = best_mae,
    coarse_landscape = coarse_landscape,
    fine_landscape   = fine_landscape,
    ambiguity        = ambiguity
  )
}


# -------------------------------------------------------------
#  recover_run_alignment
#
#  Main entry point.  Recovers recording-absolute start time(s)
#  for the paced breathing run(s) for one participant.
#
#  Arguments:
#    pipeline_trough_times  — pipelineResult$trough_times
#    recording_duration_s   — max(pipelineResult$final_time)
#    run1_trial_data        — data.frame of ses1 paced staircase trials
#                             (NULL if firstCondition == "Visual")
#    run2_trial_data        — data.frame of ses2 paced staircase trials
#                             (NULL if single-session participant)
#    run1_test_data         — data.frame of ses1 test-phase trials
#                             (from testFile.csv; NULL if unavailable).
#                             Passed through normalise_test_data() then
#                             appended to run1_trial_data in the cost
#                             function to add ~40 extra expected troughs
#                             that anchor the tail of the session.
#    run2_test_data         — same for ses2.
#    first_condition        — "Breath" or "Visual"
#    visual_session_span_s  — duration of visual ses1 task (s),
#                             i.e. max(trial.stopped) for ses1.
#                             Used as fallback lower bound when
#                             trigger_onset_ses1_s is unavailable.
#    trigger_onset_ses1_s   — recording-absolute time of the ses1
#                             start trigger (s), if available.
#                             When supplied for Visual participants,
#                             the ses2 lower bound is sharpened to
#                             trigger_onset_ses1_s +
#                             visual_session_span_s + min_gap_s
#                             instead of visual_session_span_s +
#                             min_gap_s alone, removing ambiguity
#                             from any pre-task recording time.
#    belt_quality           — "good", "degraded", or "unusable".
#                             Returns immediately with a note if
#                             "unusable" — no search is run.
#    STARTDUR               — base breath duration (default 4)
#    NUMBREATHS             — breaths per trial (default 4)
#    min_gap_s              — minimum inter-run gap (default 120)
#    coarse_step_s          — coarse search step (default 1)
#    fine_step_s            — fine search step   (default 0.1)
#    fine_window_s          — fine search half-width (default 5)
#    joint_window_s         — half-width of 2-D joint refinement
#                             for two-run participants (default 2)
#    max_match_dist_s       — match penalty cap (default 2)
#
#  Returns a list with:
#    t1, t2               — recovered ses1/ses2 start (s)
#    mae1, mae2           — MAE at best alignment per run
#    ambiguity1/2         — ambiguity score (0=clear, 1=flat)
#    n_exp1/2             — expected trough counts per run
#    n_matched1/2         — troughs matched within threshold
#    cost_landscape1/2    — coarse landscape data.frames
#    fine_landscape1/2    — fine landscape data.frames
#    notes                — character vector of warnings/flags
# -------------------------------------------------------------
recover_run_alignment <- function(pipeline_trough_times,
                                  recording_duration_s,
                                  run1_trial_data,
                                  run2_trial_data        = NULL,
                                  run1_test_data         = NULL,
                                  run2_test_data         = NULL,
                                  first_condition        = "Breath",
                                  visual_session_span_s  = NULL,
                                  trigger_onset_ses1_s   = NULL,
                                  belt_quality           = "good",
                                  STARTDUR               = 4,
                                  NUMBREATHS             = 4,
                                  min_gap_s              = 120,
                                  coarse_step_s          = 1,
                                  fine_step_s            = 0.1,
                                  fine_window_s          = 5,
                                  joint_window_s         = 2,
                                  max_match_dist_s       = 2) {

  notes <- character(0)

  # Initialise output (populated below)
  result <- list(
    t1 = NA_real_, t2 = NA_real_,
    mae1 = NA_real_, mae2 = NA_real_,
    ambiguity1 = NA_real_, ambiguity2 = NA_real_,
    n_exp1 = NA_integer_, n_exp2 = NA_integer_,
    n_matched1 = NA_integer_, n_matched2 = NA_integer_,
    cost_landscape1 = NULL, cost_landscape2 = NULL,
    fine_landscape1 = NULL, fine_landscape2 = NULL,
    notes = notes
  )

  # ------------------------------------------------------------------
  # Belt quality gate: skip search entirely for unusable signals.
  # The cost landscape over a flat/noisy belt is meaningless, and a
  # spurious "best" offset will corrupt downstream trial extraction.
  # ------------------------------------------------------------------
  if (!is.null(belt_quality) && belt_quality == "unusable") {
    result$notes <- "Alignment skipped: unusable belt signal"
    message("  [AlignRecovery] Skipping — belt quality is unusable.")
    return(result)
  }

  # ------------------------------------------------------------------
  # Helpers
  # ------------------------------------------------------------------
  run_span <- function(td) {
    if (is.null(td) || nrow(td) == 0) return(0)
    max(td$trial.stopped, na.rm = TRUE)
  }

  span1 <- run_span(run1_trial_data)
  span2 <- run_span(run2_trial_data)

  # Normalise test data to staircase column format
  test1_norm <- normalise_test_data(run1_test_data)
  test2_norm <- normalise_test_data(run2_test_data)

  # Synthesiser closures (capture trial + test data, STARTDUR, NUMBREATHS)
  synth1 <- function(t) synthesise_expected_troughs(run1_trial_data, t,
                                                     STARTDUR, NUMBREATHS,
                                                     test_data = test1_norm)
  synth2 <- function(t) synthesise_expected_troughs(run2_trial_data, t,
                                                     STARTDUR, NUMBREATHS,
                                                     test_data = test2_norm)

  n_exp1 <- if (!is.null(run1_trial_data))
    (nrow(run1_trial_data) + if (!is.null(test1_norm)) nrow(test1_norm) else 0L) *
    (NUMBREATHS + 1L) else 0L
  n_exp2 <- if (!is.null(run2_trial_data))
    (nrow(run2_trial_data) + if (!is.null(test2_norm)) nrow(test2_norm) else 0L) *
    (NUMBREATHS + 1L) else 0L

  result$n_exp1 <- n_exp1
  result$n_exp2 <- n_exp2

  count_matched <- function(synth_fn, offset, n_exp) {
    if (n_exp == 0) return(0L)
    exp_tr <- synth_fn(offset)
    sum(vapply(exp_tr, function(et)
      min(abs(pipeline_trough_times - et)) <= max_match_dist_s,
      logical(1)))
  }

  # ==================================================================
  #  CASE 1: firstCondition == "Visual"
  #  Only ses2 is paced.  Lower bound for ses2 search:
  #
  #  If trigger_onset_ses1_s is available (ses1 trigger fired):
  #    t2_min = trigger_onset_ses1_s + visual_session_span_s + min_gap_s
  #    This is recording-absolute and accounts for any pre-task delay
  #    between recording start and the first task trigger.
  #
  #  Otherwise (no ses1 trigger):
  #    t2_min = visual_session_span_s + min_gap_s
  #    This assumes recording started near task start — may be slightly
  #    too early but the search range is generous enough to compensate.
  # ==================================================================
  if (first_condition == "Visual") {

    if (is.null(run2_trial_data) || nrow(run2_trial_data) == 0) {
      result$notes <- c(notes, "Visual condition but no run2 trial data supplied.")
      return(result)
    }

    if (is.null(visual_session_span_s)) {
      notes <- c(notes,
                 "visual_session_span_s not supplied — lower bound set to min_gap_s only.")
      visual_session_span_s <- 0
    }

    # FIX: use trigger_onset_ses1_s as the recording-absolute anchor when
    # available, so the lower bound is not contaminated by pre-task
    # recording time.
    ses1_anchor <- if (!is.null(trigger_onset_ses1_s)) {
      message(sprintf(
        "  [AlignRecovery] Visual ses1 trigger at %.2f s — using as lower-bound anchor.",
        trigger_onset_ses1_s))
      trigger_onset_ses1_s
    } else {
      message("  [AlignRecovery] No ses1 trigger — lower bound uses visual_session_span_s from t=0.")
      0
    }

    t2_min <- ses1_anchor + visual_session_span_s + min_gap_s
    t2_max <- recording_duration_s - span2

    message(sprintf("  Visual/Breath: searching ses2 offset in [%.1f, %.1f] s",
                    t2_min, t2_max))

    s2 <- .two_pass_search(pipeline_trough_times, synth2,
                           t2_min, t2_max,
                           coarse_step_s, fine_step_s, fine_window_s,
                           max_match_dist_s)

    if (is.null(s2)) {
      result$notes <- c(notes, "ses2 search returned no result — range may be too narrow.")
      return(result)
    }

    result$t2              <- s2$best_offset
    result$mae2            <- s2$best_mae
    result$ambiguity2      <- s2$ambiguity
    result$cost_landscape2 <- s2$coarse_landscape
    result$fine_landscape2 <- s2$fine_landscape
    result$n_matched2      <- count_matched(synth2, s2$best_offset, n_exp2)

    if (!is.na(s2$ambiguity) && s2$ambiguity > 0.7)
      notes <- c(notes, sprintf(
        "ses2 alignment ambiguity is high (%.2f) — cost landscape may be flat.",
        s2$ambiguity))

    message(sprintf(
      "  ses2 best offset: %.2f s  |  MAE: %.3f s  |  ambiguity: %.2f  |  matched: %d/%d",
      result$t2, result$mae2, result$ambiguity2, result$n_matched2, n_exp2))

    result$notes <- notes
    return(result)
  }

  # ==================================================================
  #  CASE 2: firstCondition == "Breath"
  #  Both sessions are paced.  Sequential ses1 → ses2 search,
  #  then joint 2-D refinement.
  # ==================================================================

  if (is.null(run1_trial_data) || nrow(run1_trial_data) == 0) {
    result$notes <- c(notes, "Breath condition but no run1 trial data supplied.")
    return(result)
  }

  two_run <- !is.null(run2_trial_data) && nrow(run2_trial_data) > 0

  # ---- ses1 search ------------------------------------------------
  t1_min <- 0
  t1_max <- if (two_run)
    recording_duration_s - span1 - min_gap_s - span2
  else
    recording_duration_s - span1

  message(sprintf("  Breath ses1: searching offset in [%.1f, %.1f] s",
                  t1_min, t1_max))

  s1 <- .two_pass_search(pipeline_trough_times, synth1,
                         t1_min, t1_max,
                         coarse_step_s, fine_step_s, fine_window_s,
                         max_match_dist_s)

  if (is.null(s1)) {
    result$notes <- c(notes, "ses1 search returned no result.")
    return(result)
  }

  result$t1              <- s1$best_offset
  result$mae1            <- s1$best_mae
  result$ambiguity1      <- s1$ambiguity
  result$cost_landscape1 <- s1$coarse_landscape
  result$fine_landscape1 <- s1$fine_landscape
  result$n_matched1      <- count_matched(synth1, s1$best_offset, n_exp1)

  if (!is.na(s1$ambiguity) && s1$ambiguity > 0.7)
    notes <- c(notes, sprintf(
      "ses1 alignment ambiguity is high (%.2f) — cost landscape may be flat.",
      s1$ambiguity))

  message(sprintf(
    "  ses1 best offset: %.2f s  |  MAE: %.3f s  |  ambiguity: %.2f  |  matched: %d/%d",
    result$t1, result$mae1, result$ambiguity1, result$n_matched1, n_exp1))

  # ---- ses2 search (two-run only) ---------------------------------
  if (two_run) {

    t2_min <- s1$best_offset + span1 + min_gap_s
    t2_max <- recording_duration_s - span2

    message(sprintf("  Breath ses2: searching offset in [%.1f, %.1f] s",
                    t2_min, t2_max))

    s2 <- .two_pass_search(pipeline_trough_times, synth2,
                           t2_min, t2_max,
                           coarse_step_s, fine_step_s, fine_window_s,
                           max_match_dist_s)

    if (is.null(s2)) {
      notes <- c(notes, "ses2 search returned no result.")
    } else {

      result$t2              <- s2$best_offset
      result$mae2            <- s2$best_mae
      result$ambiguity2      <- s2$ambiguity
      result$cost_landscape2 <- s2$coarse_landscape
      result$fine_landscape2 <- s2$fine_landscape
      result$n_matched2      <- count_matched(synth2, s2$best_offset, n_exp2)

      if (!is.na(s2$ambiguity) && s2$ambiguity > 0.7)
        notes <- c(notes, sprintf(
          "ses2 alignment ambiguity is high (%.2f) — cost landscape may be flat.",
          s2$ambiguity))

      message(sprintf(
        "  ses2 best offset: %.2f s  |  MAE: %.3f s  |  ambiguity: %.2f  |  matched: %d/%d",
        result$t2, result$mae2, result$ambiguity2, result$n_matched2, n_exp2))

      # ---- Joint 2-D refinement -----------------------------------
      # Fine grid search over a small window around both sequential
      # best offsets, enforcing the t2 > t1 + span1 + min_gap constraint.
      # Weighted combined MAE = (mae1*n_exp1 + mae2*n_exp2) / (n_exp1+n_exp2).
      message("  Running joint 2-D refinement...")

      t1_range <- seq(max(t1_min, result$t1 - joint_window_s),
                      min(t1_max, result$t1 + joint_window_s),
                      by = fine_step_s)
      t2_range <- seq(max(t2_min, result$t2 - joint_window_s),
                      min(t2_max, result$t2 + joint_window_s),
                      by = fine_step_s)

      # Pre-compute ses1 MAEs across t1_range to avoid redundant calls
      # in the inner loop — each t1 candidate is independent of t2.
      mae1_vec <- vapply(t1_range, function(t1c)
        score_trough_alignment(pipeline_trough_times, synth1(t1c),
                               max_match_dist_s),
        numeric(1))

      best_joint_mae <- Inf
      best_t1_joint  <- result$t1
      best_t2_joint  <- result$t2

      for (i_t1 in seq_along(t1_range)) {
        t1c    <- t1_range[i_t1]
        mae1c  <- mae1_vec[i_t1]
        for (t2c in t2_range) {
          if (t2c < t1c + span1 + min_gap_s) next
          mae2c     <- score_trough_alignment(pipeline_trough_times,
                                              synth2(t2c), max_match_dist_s)
          joint_mae <- (mae1c * n_exp1 + mae2c * n_exp2) / (n_exp1 + n_exp2)
          if (joint_mae < best_joint_mae) {
            best_joint_mae <- joint_mae
            best_t1_joint  <- t1c
            best_t2_joint  <- t2c
          }
        }
      }

      prev_joint_mae <- (result$mae1 * n_exp1 + result$mae2 * n_exp2) /
        (n_exp1 + n_exp2)

      if (best_joint_mae < prev_joint_mae) {
        message(sprintf("  Joint refinement improved combined MAE: %.3f → %.3f",
                        prev_joint_mae, best_joint_mae))
        result$t1         <- best_t1_joint
        result$t2         <- best_t2_joint
        result$mae1       <- score_trough_alignment(pipeline_trough_times,
                                                    synth1(best_t1_joint),
                                                    max_match_dist_s)
        result$mae2       <- score_trough_alignment(pipeline_trough_times,
                                                    synth2(best_t2_joint),
                                                    max_match_dist_s)
        result$n_matched1 <- count_matched(synth1, best_t1_joint, n_exp1)
        result$n_matched2 <- count_matched(synth2, best_t2_joint, n_exp2)
      } else {
        message("  Joint refinement: no improvement over sequential result.")
      }
    }
  }

  # ------------------------------------------------------------------
  # POST-SEARCH GUARD: ses2 onset must not fall inside the ses1 window
  # ------------------------------------------------------------------
  # Failure mode: when the ses2 signal is entirely unusable (e.g. belt
  # flatline, participant removed belt between sessions), the ses2 cost
  # function finds no real troughs and latches onto the ses1 trough
  # cluster as the best match — returning a ses2 onset nearly identical
  # to ses1.  This produces a spurious alignment that is undetectable
  # from the MAE or ambiguity score alone.
  #
  # Guard: if the recovered ses2 onset falls within the ses1 task window
  # (t2 < t1 + run1_span + min_gap_s / 2), the onset is almost certainly
  # spurious.  Set t2 to NA and add a warning note so the main loop knows
  # to treat ses2 physio as unavailable.
  #
  # The half-gap tolerance (min_gap_s / 2 = 60 s) gives a buffer for
  # genuine short inter-session gaps while still catching the pathological
  # case where t2 ≈ t1 (difference < 1 s, as seen in P13081).
  # ------------------------------------------------------------------
  if (!is.na(result$t2) && !is.null(run1_trial_data)) {
    ses1_span_for_guard <- max(run1_trial_data$trial.stopped, na.rm = TRUE)

    # Anchor on the best available ses1 onset — trigger takes priority over
    # recovery, because recovery can be wrong (e.g. P13081: recovery found
    # 157s while trigger correctly says 2211s). Using a wrong recovery t1
    # as the anchor would make the guard threshold far too low, allowing a
    # spurious ses2 onset at t1+~2000s to pass unchecked.
    ses1_anchor_for_guard <- if (!is.null(trigger_onset_ses1_s) &&
                                  !is.na(trigger_onset_ses1_s))
      trigger_onset_ses1_s
    else if (!is.na(result$t1))
      result$t1
    else
      NA_real_

    guard_threshold <- if (!is.na(ses1_anchor_for_guard))
      ses1_anchor_for_guard + ses1_span_for_guard + min_gap_s / 2
    else
      NA_real_

    if (!is.na(guard_threshold) && result$t2 < guard_threshold) {
      msg <- sprintf(
        paste0("ses2 onset (%.1f s) falls within ses1 window ",
               "(t1=%.1f s + span=%.1f s + half-gap=%.1f s = %.1f s). ",
               "Likely caused by flatlined/absent ses2 signal latching onto ",
               "ses1 troughs. ses2 onset set to NA."),
        result$t2, result$t1, ses1_span_for_guard,
        min_gap_s / 2, guard_threshold)
      message(sprintf("  [WARN] %s", msg))
      notes      <- c(notes, msg)
      result$t2         <- NA_real_
      result$mae2       <- NA_real_
      result$ambiguity2 <- NA_real_
      result$n_matched2 <- NA_integer_
    }
  }

  result$notes <- notes
  result
}


# -------------------------------------------------------------
#  plot_alignment_cost
#
#  Diagnostic plot of the MAE cost landscape(s).
#  Essential for manual review — a well-identified alignment
#  shows a clear single trough; ambiguous cases show multiple
#  local minima or a flat landscape.
#
#  Arguments:
#    recovery_result — output of recover_run_alignment()
#    participant_id  — label for plot title
#    plot_file       — path to save PNG, or NULL to display
# -------------------------------------------------------------
plot_alignment_cost <- function(recovery_result,
                                participant_id = "",
                                plot_file      = NULL) {
  library(ggplot2)
  library(patchwork)

  make_panel <- function(landscape, fine_landscape, best_t, mae_val,
                         amb_val, run_label) {
    if (is.null(landscape)) return(NULL)

    p_coarse <- ggplot(landscape, aes(x = offset, y = mae)) +
      geom_line(colour = "#2c7be5", linewidth = 0.6) +
      geom_vline(xintercept = best_t,
                 colour = "darkorange", linetype = "dashed", linewidth = 0.8) +
      annotate("text", x = best_t, y = max(landscape$mae, na.rm = TRUE),
               label = sprintf("t = %.2f s\nMAE = %.3f s\nambiguity = %.2f",
                               best_t, mae_val, amb_val),
               hjust = -0.1, vjust = 1, size = 3, colour = "darkorange") +
      labs(title = sprintf("%s — coarse landscape", run_label),
           x = "Candidate offset (s)", y = "MAE (s)") +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(size = 9, face = "bold"))

    p_fine <- ggplot(fine_landscape, aes(x = offset, y = mae)) +
      geom_line(colour = "#e63946", linewidth = 0.6) +
      geom_vline(xintercept = best_t,
                 colour = "darkorange", linetype = "dashed", linewidth = 0.8) +
      labs(title = sprintf("%s — fine landscape (±window around best)", run_label),
           x = "Candidate offset (s)", y = "MAE (s)") +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(size = 9, face = "bold"))

    list(coarse = p_coarse, fine = p_fine)
  }

  r  <- recovery_result
  p1 <- make_panel(r$cost_landscape1, r$fine_landscape1,
                   r$t1, r$mae1, r$ambiguity1, "Session 1")
  p2 <- make_panel(r$cost_landscape2, r$fine_landscape2,
                   r$t2, r$mae2, r$ambiguity2, "Session 2")

  panels <- c(
    if (!is.null(p1)) list(p1$coarse, p1$fine) else NULL,
    if (!is.null(p2)) list(p2$coarse, p2$fine) else NULL
  )

  if (length(panels) == 0) {
    message("No landscapes to plot.")
    return(invisible(NULL))
  }

  title_str <- paste0(
    "Alignment Recovery",
    if (nchar(participant_id) > 0) paste0("  —  P", participant_id) else "",
    if (length(r$notes) > 0)
      paste0("\n\u26a0 ", paste(r$notes, collapse = "\n\u26a0 ")) else ""
  )

  combined <- patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(
      title = title_str,
      theme = theme(plot.title = element_text(
        size   = 11, face = "bold",
        colour = if (length(r$notes) > 0) "#e63946" else "black")))

  if (!is.null(plot_file)) {
    n_rows <- ceiling(length(panels) / 2)
    ggsave(plot_file, combined, width = 14, height = 4 * n_rows, dpi = 120)
    message("Alignment plot saved to: ", plot_file)
  } else {
    print(combined)
  }

  invisible(combined)
}


# -------------------------------------------------------------
#  plot_alignment_overview
#
#  Four-row per-participant QC figure combining:
#    Row 1 (full-width): Full recording signal with detected
#            troughs and shaded session windows.
#    Row 2 (two panels): Cost landscapes (coarse + fine) for
#            each session. Visual ses1 shows a text note.
#    Row 3 (two panels): Detected trough count per trial,
#            against the expected 5. Visual ses1 shows a note.
#    Row 4 (two panels): Signal zoom over the first
#            n_zoom_trials trials, with observed troughs
#            (green triangles) and expected trough positions
#            (coloured tick marks below signal).
#
#  Visual-first participants: Session 1 panels display a note
#  explaining that pacing was not expected and no alignment was
#  attempted.
#
#  Arguments:
#    processedSignal  — pipelineResult$final_signal
#    processedTime    — pipelineResult$final_time
#    trough_times     — pipelineResult$trough_times
#    recovery         — output of recover_run_alignment()
#    run1_trial_data  — ses1 staircase trials (NULL if Visual)
#    run2_trial_data  — ses2 staircase trials
#    run1_test_data   — ses1 test trials, already normalised via
#                       normalise_test_data() (NULL if unavailable)
#    run2_test_data   — ses2 test trials, already normalised
#    runOnsetSec      — numeric(2): recording-absolute onsets;
#                       NA where unavailable
#    first_condition  — "Breath" or "Visual"
#    participant_id   — label for plot titles
#    LAG              — hardware latency correction (s)
#    STARTDUR / NUMBREATHS — passed to synthesise_expected_troughs
#    n_zoom_trials    — trials shown in the zoom panel (default 6)
#    save_path        — file path for PNG, or NULL to print
#    width / height   — plot dimensions in inches
# -------------------------------------------------------------
plot_alignment_overview <- function(processedSignal,
                                    processedTime,
                                    trough_times,
                                    recovery         = NULL,
                                    run1_trial_data  = NULL,
                                    run2_trial_data  = NULL,
                                    run1_test_data   = NULL,
                                    run2_test_data   = NULL,
                                    runOnsetSec      = c(NA_real_, NA_real_),
                                    first_condition  = "Breath",
                                    participant_id   = "",
                                    LAG              = 0.2,
                                    STARTDUR         = 4,
                                    NUMBREATHS       = 4,
                                    n_zoom_trials    = 6,
                                    save_path        = NULL,
                                    width            = 16,
                                    height           = 14) {
  library(ggplot2)
  library(patchwork)

  C_BLUE   <- "#2c7be5"; C_RED <- "#e63946"
  C_GREEN  <- "#2dc653"; C_ORANGE <- "darkorange"; C_GREY <- "#adb5bd"
  ses_cols <- c("1" = C_ORANGE, "2" = C_RED)

  # ── Helper: count troughs in padded trial window ──────────
  count_trial_troughs <- function(onset, trial_df,
                                  start_col = "trial.started",
                                  stop_col  = "trial.stopped",
                                  pad = 1.0) {
    vapply(seq_len(nrow(trial_df)), function(i) {
      ss <- onset + LAG + trial_df[[start_col]][i]
      se <- onset + LAG + trial_df[[stop_col]][i]
      sum(trough_times >= ss - pad & trough_times <= se + pad)
    }, integer(1))
  }

  # ── Helper: placeholder panel with centred text ───────────
  blank_panel <- function(title, msg, col = "grey50") {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = msg,
               size = 3.2, colour = col, hjust = 0.5, vjust = 0.5) +
      theme_void() +
      labs(title = title) +
      theme(plot.title = element_text(size = 9, face = "bold"))
  }

  # ── ROW 1: Full recording ──────────────────────────────────
  ds_step <- max(1L, floor(length(processedSignal) / 8000L))
  idx_ds  <- seq(1L, length(processedTime), ds_step)
  df_full <- data.frame(time   = processedTime[idx_ds],
                         signal = processedSignal[idx_ds])

  tr_samp  <- pmax(1L, pmin(length(processedSignal),
                             round(trough_times *
                                     1 / (processedTime[2] - processedTime[1])) + 1L))
  df_tr    <- data.frame(time   = trough_times,
                          signal = processedSignal[tr_samp])

  p_full <- ggplot(df_full, aes(time, signal)) +
    geom_line(colour = C_BLUE, linewidth = 0.3, alpha = 0.7) +
    geom_point(data = df_tr, aes(time, signal),
               shape = 25, fill = C_GREEN, colour = C_GREEN,
               size = 1.0, alpha = 0.4)

  for (run_i in 1:2) {
    onset   <- runOnsetSec[run_i]
    trial_d <- if (run_i == 1) run1_trial_data else run2_trial_data
    if (!is.na(onset) && !is.null(trial_d) && nrow(trial_d) > 0) {
      span <- max(trial_d$trial.stopped, na.rm = TRUE)
      col  <- ses_cols[as.character(run_i)]
      is_vis1 <- (run_i == 1 && first_condition == "Visual")
      sig_top  <- quantile(processedSignal, 0.98, na.rm = TRUE)
      p_full <- p_full +
        annotate("rect", xmin = onset, xmax = onset + span,
                 ymin = -Inf, ymax = Inf, fill = col, alpha = 0.07) +
        annotate("segment", x = onset, xend = onset, y = -Inf, yend = Inf,
                 colour = col, linewidth = 0.9,
                 linetype = if (is_vis1) "dotted" else "dashed") +
        annotate("label", x = onset + span * 0.02, y = sig_top,
                 label = if (is_vis1) "Session 1\n(Visual — no pacing)"
                         else sprintf("Session %d\n(Breath)", run_i),
                 colour = col, size = 2.6, fontface = "bold",
                 label.size = 0.15, hjust = 0, vjust = 1)
    }
  }

  p_full <- p_full +
    labs(title = sprintf("P%s — full recording with recovered onsets", participant_id),
         x = "Time (s)", y = "z-score") +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(size = 10, face = "bold"))

  # ── ROW 2: Cost landscapes ─────────────────────────────────
  make_landscape <- function(run_i) {
    is_vis1 <- (run_i == 1 && first_condition == "Visual")
    if (is_vis1)
      return(blank_panel(sprintf("Session %d — cost landscape", run_i),
                         "Session 1 is Visual\nalignment not performed"))
    if (is.null(recovery)) return(blank_panel(
      sprintf("Session %d — cost landscape", run_i), "No alignment data"))

    coarse <- if (run_i==1) recovery$cost_landscape1 else recovery$cost_landscape2
    fine   <- if (run_i==1) recovery$fine_landscape1 else recovery$fine_landscape2
    best   <- if (run_i==1) recovery$t1    else recovery$t2
    mae    <- if (run_i==1) recovery$mae1  else recovery$mae2
    amb    <- if (run_i==1) recovery$ambiguity1 else recovery$ambiguity2
    col    <- ses_cols[as.character(run_i)]

    if (is.null(coarse) || is.na(best))
      return(blank_panel(sprintf("Session %d — cost landscape", run_i),
                         "No landscape data"))

    ymax <- max(fine$mae, na.rm = TRUE)
    ggplot() +
      geom_line(data = coarse, aes(offset, mae), colour = C_GREY, linewidth = 0.6) +
      geom_line(data = fine,   aes(offset, mae), colour = col,    linewidth = 1.0) +
      geom_vline(xintercept = best, linetype = "dashed",
                 colour = "black", linewidth = 0.8) +
      annotate("label", x = best, y = ymax,
               label = sprintf("t = %.1f s\nMAE = %.3f s\namb = %.2f",
                                best, mae, amb),
               size = 2.5, hjust = -0.1, vjust = 1,
               label.size = 0.15, colour = "black") +
      labs(title = sprintf("Session %d — cost landscape", run_i),
           x = "Candidate onset (s)", y = "MAE (s)") +
      theme_minimal(base_size = 9) +
      theme(plot.title = element_text(size = 9, face = "bold"))
  }

  # ── ROW 3: Trough counts ───────────────────────────────────
  make_count <- function(run_i, trial_df) {
    onset   <- runOnsetSec[run_i]
    col     <- ses_cols[as.character(run_i)]
    is_vis1 <- (run_i == 1 && first_condition == "Visual")

    if (is_vis1)
      return(blank_panel(sprintf("Session %d — trough count per trial", run_i),
                         "Session 1 is Visual\ntrough count not meaningful"))
    if (is.null(trial_df) || nrow(trial_df) == 0 || is.na(onset))
      return(blank_panel(sprintf("Session %d — trough count per trial", run_i),
                         "No trial / onset data"))

    counts <- pmin(count_trial_troughs(onset, trial_df), 8L)
    df_c   <- data.frame(trial = seq_along(counts), n = counts)
    exp_n  <- NUMBREATHS + 1L

    ggplot(df_c, aes(trial, n)) +
      geom_col(fill = col, alpha = 0.7, width = 0.75) +
      geom_hline(yintercept = exp_n, linetype = "dashed",
                 colour = "black", linewidth = 0.7) +
      annotate("text", x = nrow(df_c) - 0.5, y = exp_n + 0.4,
               label = sprintf("Expected (%d)", exp_n),
               size = 2.4, hjust = 1, colour = "black") +
      scale_y_continuous(limits = c(0, 8), breaks = 0:8) +
      labs(title = sprintf("Session %d — detected troughs per trial (±1 s)", run_i),
           x = "Trial", y = "Troughs") +
      theme_minimal(base_size = 9) +
      theme(plot.title = element_text(size = 9, face = "bold"))
  }

  # ── ROW 4: Signal zoom ─────────────────────────────────────
  make_zoom <- function(run_i, trial_df, test_df) {
    onset   <- runOnsetSec[run_i]
    col     <- ses_cols[as.character(run_i)]
    is_vis1 <- (run_i == 1 && first_condition == "Visual")
    dt      <- processedTime[2] - processedTime[1]

    if (is.null(trial_df) || nrow(trial_df) == 0 || is.na(onset))
      return(blank_panel(sprintf("Session %d — trial zoom", run_i),
                         "No trial / onset data"))

    zoom_df <- head(trial_df, n_zoom_trials)
    t_lo <- onset + LAG + zoom_df$trial.started[1] - 5
    t_hi <- onset + LAG + zoom_df$trial.stopped[nrow(zoom_df)] + 5

    # Signal segment
    sig_mask <- processedTime >= t_lo & processedTime <= t_hi
    df_sig   <- data.frame(time   = processedTime[sig_mask],
                            signal = processedSignal[sig_mask])

    # Observed troughs in window
    tr_mask <- trough_times >= t_lo & trough_times <= t_hi
    tr_samp_z <- pmax(1L, pmin(length(processedSignal),
                                round(trough_times[tr_mask] / dt) + 1L))
    df_tr_z <- data.frame(time   = trough_times[tr_mask],
                           signal = processedSignal[tr_samp_z])

    # Expected troughs
    exp_all  <- synthesise_expected_troughs(trial_df, onset,
                                            STARTDUR, NUMBREATHS, test_df)
    exp_mask <- exp_all >= t_lo & exp_all <= t_hi
    df_exp   <- data.frame(time = exp_all[exp_mask])

    # Trial rectangles
    df_rect <- data.frame(xmin = onset + LAG + zoom_df$trial.started,
                           xmax = onset + LAG + zoom_df$trial.stopped)

    y_lo <- min(df_sig$signal, na.rm = TRUE)
    tick_lo <- y_lo - 0.3
    tick_hi <- y_lo

    p <- ggplot() +
      geom_rect(data = df_rect,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                fill = col, alpha = 0.07, inherit.aes = FALSE) +
      geom_vline(data = df_rect, aes(xintercept = xmin),
                 colour = col, linewidth = 0.4, alpha = 0.5) +
      geom_line(data = df_sig, aes(time, signal),
                colour = C_BLUE, linewidth = 0.7, alpha = 0.85)

    if (nrow(df_tr_z) > 0)
      p <- p + geom_point(data = df_tr_z, aes(time, signal),
                           shape = 25, fill = C_GREEN, colour = C_GREEN,
                           size = 3.5, zorder = 4)

    if (nrow(df_exp) > 0)
      p <- p + geom_segment(data = df_exp,
                              aes(x = time, xend = time,
                                  y = tick_lo, yend = tick_hi),
                              colour = col, linewidth = 1.1, alpha = 0.9)

    label_suffix <- if (is_vis1) sprintf(
      "Session %d (Visual) — first %d trials\n[no pacing expected]",
      run_i, n_zoom_trials)
    else
      sprintf("Session %d — first %d trial windows  ▼obs  %s exp",
              run_i, n_zoom_trials,
              if (run_i==1) "▬" else "▬")

    p + labs(title = label_suffix, x = "Time (s)", y = "z-score") +
      theme_minimal(base_size = 9) +
      theme(plot.title = element_text(size = 9, face = "bold"))
  }

  # ── Assemble ─────────────────────────────────────────────────
  eff_fs <- if (length(processedTime) > 1)
    round(1 / (processedTime[2] - processedTime[1]), 2) else NA

  combined <- p_full /
    (make_landscape(1) | make_landscape(2)) /
    (make_count(1, run1_trial_data) | make_count(2, run2_trial_data)) /
    (make_zoom(1, run1_trial_data, run1_test_data) |
     make_zoom(2, run2_trial_data, run2_test_data)) +
    patchwork::plot_layout(heights = c(1.8, 1, 1, 1.4)) +
    patchwork::plot_annotation(
      title = sprintf(
        "Alignment Overview — P%s  |  ses1: %s  |  fs: %.2f Hz",
        participant_id, first_condition, eff_fs),
      theme = theme(plot.title = element_text(size = 11, face = "bold")))

  if (!is.null(save_path)) {
    ggsave(save_path, combined, width = width, height = height, dpi = 120)
    message(sprintf("  [AlignOverview] Saved: %s", basename(save_path)))
  } else {
    print(combined)
  }

  invisible(combined)
}
