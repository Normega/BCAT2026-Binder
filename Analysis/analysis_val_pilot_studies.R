# ============================================================
# analysis_val_pilot_studies.R
# BCAT — Task Validation: Pilot Study Comparisons and Convergence
#
# Corresponds to: Supplement S1.1 (Study 1A/1B), S1.2 (convergence),
#                 S1.3 (Study 3 salience-magnitude orthogonality)
#
# Assumes in environment (set by MainAnalysis.R):
#   s1l, s1s, s2l, s3l
#   RESULTS_DIR
#
# Outputs:
#   study1_task_comparison.csv          — S1.1: MRC stability
#   table_study3_salience_accuracy.csv  — S1.3: salience orthogonality
#   table_staircase_convergence.csv     — S1.2 + S1.4: convergence cells
#   study2_convergence_elbow.csv        — S1.2: per-trial cumulative SD
#   study2_convergence_elbow_summary.csv — S1.2: elbow summary
# ============================================================

message("\n========================================")
message("TASK VALIDATION: Pilot studies, convergence, salience orthogonality")
message("========================================")


# ── S1.1 Study 1 task comparison: MRC stability ──────────────
#
# Both tasks produce comparable MRC estimates, confirmed by
# cross-task correlation. The adaptive staircase (Task A) is
# substantially more stable trial-to-trial than the ascending
# limits approach (Task B), motivating adoption of Quest for
# Studies 4 and 5.
#
# MRC definitions:
#   Task A: mean of last 4 TotalRateChange values per participant
#   Task B: mean abs(DetectedChange - 1) on accurate change trials

mrc_A <- s1l_A_sorted |>
  dplyr::group_by(id) |>
  dplyr::summarise(
    MRC_A = mean(utils::tail(TotalRateChange, 4)),
    SD_A  = stats::sd(utils::tail(TotalRateChange, 4)),
    .groups = "drop"
  )

mrc_B <- s1l |>
  dplyr::filter(
    Group     == "TaskB",
    Direction %in% c("Faster", "Slower"),
    Accuracy  == 1
  ) |>
  dplyr::mutate(thresh_deviation = abs(as.numeric(DetectedChange) - 1)) |>
  dplyr::group_by(id) |>
  dplyr::summarise(
    MRC_B        = mean(thresh_deviation, na.rm = TRUE),
    SD_B         = stats::sd(thresh_deviation,  na.rm = TRUE),
    n_acc_trials = dplyr::n(),
    .groups      = "drop"
  )

mrc_compare <- dplyr::inner_join(mrc_A, mrc_B, by = "id") |>
  tidyr::drop_na()

mrc_cor      <- stats::cor.test(mrc_compare$MRC_A, mrc_compare$MRC_B)
cv_A         <- mean(mrc_compare$SD_A, na.rm = TRUE) /
                mean(mrc_compare$MRC_A, na.rm = TRUE)
cv_B         <- mean(mrc_compare$SD_B, na.rm = TRUE) /
                mean(mrc_compare$MRC_B, na.rm = TRUE)
stability_t  <- stats::t.test(mrc_compare$SD_A, mrc_compare$SD_B, paired = TRUE)

cat("\n--- Study 1 task comparison: MRC and stability ---\n")
cat(sprintf("MRC_A: M = %.3f  MRC_B: M = %.3f\n",
            mean(mrc_compare$MRC_A), mean(mrc_compare$MRC_B)))
cat(sprintf("Cross-task correlation: r = %.3f, p = %.4f, N = %d\n",
            mrc_cor$estimate, mrc_cor$p.value, nrow(mrc_compare)))
cat(sprintf("CV Task A: %.3f  CV Task B: %.3f\n", cv_A, cv_B))
cat(sprintf("Paired t on SD: t(%d) = %.3f, p = %.4f\n",
            stability_t$parameter, stability_t$statistic,
            stability_t$p.value))

readr::write_csv(
  tibble::tibble(
    measure = c("MRC_A_mean", "MRC_B_mean",
                "cross_task_r", "cross_task_p",
                "CV_TaskA", "CV_TaskB",
                "stability_t", "stability_p", "N"),
    value   = c(mean(mrc_compare$MRC_A),
                mean(mrc_compare$MRC_B),
                mrc_cor$estimate,
                mrc_cor$p.value,
                cv_A, cv_B,
                stability_t$statistic,
                stability_t$p.value,
                nrow(mrc_compare))
  ),
  file.path(RESULTS_DIR, "study1_task_comparison.csv")
)
message("Saved: study1_task_comparison.csv")


# ── S1.3 Study 3: salience effect on accuracy ────────────────
#
# Confirms salience and magnitude are orthogonal: identical
# magnitude distributions in both salience conditions, yet
# high-salience changes detected at more than twice the rate.

s3_sal_acc <- s3l |>
  dplyr::filter(!is.na(Accuracy), !is.na(Salience)) |>
  dplyr::group_by(id, Salience) |>
  dplyr::summarise(acc = mean(Accuracy), .groups = "drop") |>
  tidyr::pivot_wider(names_from = Salience, values_from = acc) |>
  tidyr::drop_na()

s3_sal_t <- stats::t.test(s3_sal_acc$High, s3_sal_acc$Low, paired = TRUE)
s3_sal_d <- (mean(s3_sal_acc$High) - mean(s3_sal_acc$Low)) /
            sd(c(s3_sal_acc$High, s3_sal_acc$Low))

readr::write_csv(
  tibble::tibble(
    study    = "Study3",
    n_pp     = nrow(s3_sal_acc),
    M_High   = mean(s3_sal_acc$High),
    M_Low    = mean(s3_sal_acc$Low),
    t        = unname(s3_sal_t$statistic),
    df       = unname(s3_sal_t$parameter),
    p        = s3_sal_t$p.value,
    cohens_d = s3_sal_d
  ),
  file.path(RESULTS_DIR, "table_study3_salience_accuracy.csv")
)
message("Saved: table_study3_salience_accuracy.csv")


# ── S1.2 + S1.4 Staircase convergence (Studies 2, 4, 5) ──────
#
# Convergence = within-condition SD of Change estimates,
# first 6 vs. last 6 trials. For Studies 4 and 5, assessed
# within each of the 4 Salience × Direction cells independently.
# Cohen's d_z = t / sqrt(n).

.convergence_cells <- function(long_data, study_label,
                                change_col   = "Change",
                                group_filter = NULL,
                                group_col    = NULL,
                                ses_col      = NULL) {
  d <- long_data
  if (!is.null(group_filter) && !is.null(group_col))
    d <- dplyr::filter(d, .data[[group_col]] == group_filter)

  d <- dplyr::filter(d, Direction %in% c("Faster", "Slower"))

  group_cols <- c("id",
                  if (!is.null(ses_col) && ses_col %in% names(d)) ses_col,
                  if ("Salience"  %in% names(d)) "Salience",
                  if ("Direction" %in% names(d)) "Direction")

  conv <- d |>
    dplyr::arrange(dplyr::across(dplyr::all_of(c(group_cols, "Trial")))) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n_trials  = dplyr::n(),
      SD_first6 = stats::sd(utils::head(.data[[change_col]], 6), na.rm = TRUE),
      SD_last6  = stats::sd(utils::tail(.data[[change_col]], 6), na.rm = TRUE),
      .groups   = "drop"
    ) |>
    dplyr::filter(n_trials >= 6) |>
    tidyr::drop_na(SD_first6, SD_last6)

  cells <- conv |>
    dplyr::select(dplyr::any_of(c("Salience", "Direction"))) |>
    dplyr::distinct() |>
    tidyr::unite("cell", dplyr::everything(), sep = "x", remove = FALSE)

  purrr::pmap_dfr(cells, function(cell, Salience = NULL, Direction = NULL) {
    sub <- conv
    if (!is.null(Salience)  && "Salience"  %in% names(conv))
      sub <- dplyr::filter(sub, Salience  == !!Salience)
    if (!is.null(Direction) && "Direction" %in% names(conv))
      sub <- dplyr::filter(sub, Direction == !!Direction)
    if (nrow(sub) < 2) return(NULL)

    tt  <- stats::t.test(sub$SD_first6, sub$SD_last6, paired = TRUE)
    dz  <- unname(tt$statistic) / sqrt(nrow(sub))
    pct <- (mean(sub$SD_first6) - mean(sub$SD_last6)) /
            mean(sub$SD_first6) * 100

    tibble::tibble(
      study         = study_label,
      cell          = cell,
      n             = nrow(sub),
      n_trials_mean = mean(sub$n_trials),
      SD_first6     = mean(sub$SD_first6),
      SD_last6      = mean(sub$SD_last6),
      pct_reduction = pct,
      t             = unname(tt$statistic),
      df            = unname(tt$parameter),
      p             = tt$p.value,
      cohens_dz     = dz
    )
  })
}

# Study 2: single staircase, no Salience/Direction breakdown
s2_conv_overall <- s2l |>
  dplyr::filter(Direction %in% c("Faster", "Slower")) |>
  dplyr::arrange(id, TrialNum) |>
  dplyr::group_by(id) |>
  dplyr::summarise(
    n_trials  = dplyr::n(),
    SD_first6 = stats::sd(utils::head(TotalRateChange, 6), na.rm = TRUE),
    SD_last6  = stats::sd(utils::tail(TotalRateChange, 6), na.rm = TRUE),
    .groups   = "drop"
  ) |> tidyr::drop_na()

s2_tt  <- stats::t.test(s2_conv_overall$SD_first6,
                         s2_conv_overall$SD_last6, paired = TRUE)
s2_row <- tibble::tibble(
  study         = "Study2", cell = "Overall",
  n             = nrow(s2_conv_overall),
  n_trials_mean = mean(s2_conv_overall$n_trials),
  SD_first6     = mean(s2_conv_overall$SD_first6),
  SD_last6      = mean(s2_conv_overall$SD_last6),
  pct_reduction = (mean(s2_conv_overall$SD_first6) -
                   mean(s2_conv_overall$SD_last6)) /
                   mean(s2_conv_overall$SD_first6) * 100,
  t             = unname(s2_tt$statistic),
  df            = unname(s2_tt$parameter),
  p             = s2_tt$p.value,
  cohens_dz     = unname(s2_tt$statistic) / sqrt(nrow(s2_conv_overall))
)

s4_conv <- .convergence_cells(
  s4l, "Study4", change_col = "Change",
  group_filter = "Breath", group_col = "Group")

s5_conv <- .convergence_cells(
  s5l, "Study5", change_col = "Change",
  group_filter = "breath", group_col = "Condition",
  ses_col = "ses")

table_staircase_convergence <- dplyr::bind_rows(s2_row, s4_conv, s5_conv) |>
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4)))

readr::write_csv(table_staircase_convergence,
                 file.path(RESULTS_DIR, "table_staircase_convergence.csv"))
message("Saved: table_staircase_convergence.csv")

s45_dz  <- table_staircase_convergence |>
  dplyr::filter(study != "Study2") |>
  dplyr::pull(cohens_dz)
s45_pct <- table_staircase_convergence |>
  dplyr::filter(study != "Study2") |>
  dplyr::pull(pct_reduction)
cat(sprintf(
  "\nConvergence summary (Studies 4/5): d_z = %.2f-%.2f; reduction = %.1f-%.1f%%\n",
  min(s45_dz), max(s45_dz), min(s45_pct), max(s45_pct)))
cat(sprintf(
  "Study 2: d_z = %.2f; reduction = %.1f%%\n",
  s2_row$cohens_dz, s2_row$pct_reduction))


# ── S1.2 Study 2 elbow analysis ──────────────────────────────
#
# Cumulative within-person SD per trial; elbow = first trial
# where rate of SD reduction drops below 10% of initial rate.

s2_change_sorted <- s2l |>
  dplyr::filter(Direction %in% c("Faster", "Slower"),
                !is.na(TotalRateChange)) |>
  dplyr::arrange(id, TrialNum)

max_t <- s2_change_sorted |>
  dplyr::group_by(id) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
  dplyr::summarise(median_n = stats::median(n)) |>
  dplyr::pull(median_n) |>
  as.integer()

running_within_sd <- purrr::map_dfr(2:max_t, function(t) {
  sds <- s2_change_sorted |>
    dplyr::group_by(id) |>
    dplyr::filter(dplyr::row_number() <= t) |>
    dplyr::summarise(
      sd_t = stats::sd(TotalRateChange, na.rm = TRUE),
      n    = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::filter(n >= 2, !is.na(sd_t))
  tibble::tibble(trial = t, mean_sd = mean(sds$sd_t), n_pp = nrow(sds))
})

running_within_sd <- running_within_sd |>
  dplyr::mutate(delta = c(NA_real_, diff(mean_sd)))

initial_rate    <- abs(running_within_sd$delta[running_within_sd$trial == 3])
elbow_threshold <- 0.10 * initial_rate

elbow_trial <- running_within_sd |>
  dplyr::filter(!is.na(delta), abs(delta) < elbow_threshold) |>
  dplyr::slice(1) |>
  dplyr::pull(trial)

cat(sprintf("\nStudy 2 convergence elbow: trial %d\n", elbow_trial))
cat(sprintf("(SD change rate drops below %.1f%% of initial rate)\n", 10))
cat(sprintf("SD at elbow: %.4f; SD at final trial: %.4f\n",
            running_within_sd$mean_sd[running_within_sd$trial == elbow_trial],
            running_within_sd$mean_sd[nrow(running_within_sd)]))

readr::write_csv(running_within_sd,
                 file.path(RESULTS_DIR, "study2_convergence_elbow.csv"))
message("Saved: study2_convergence_elbow.csv")

readr::write_csv(
  tibble::tibble(
    study               = "Study2",
    elbow_trial         = elbow_trial,
    elbow_threshold_pct = 10,
    sd_at_elbow         = running_within_sd$mean_sd[
      running_within_sd$trial == elbow_trial],
    sd_at_final         = running_within_sd$mean_sd[nrow(running_within_sd)],
    n_change_trials_median = max_t,
    note = paste0("Elbow = first trial where |delta SD| < 10% of initial rate. ",
                  "~17 change trials from ~25 total (1/3 no-change interleaved). ",
                  "Matches ~8 effective change trials per condition in Studies 4/5.")
  ),
  file.path(RESULTS_DIR, "study2_convergence_elbow_summary.csv")
)
message("Saved: study2_convergence_elbow_summary.csv")
