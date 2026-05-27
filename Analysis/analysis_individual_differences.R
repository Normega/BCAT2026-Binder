# ============================================================
# analysis_individual_differences.R
# BCAT — Individual Difference Feature Extraction
#
# Extracts participant-level features combining task performance
# and within-person relationships for correlation with questionnaire
# measures. Primary focus: Studies 4 and 5 (replication design).
# Studies 1A, 2, 3 included where features are computable.
#
# Output:
#   id_features_study{N}.csv  — per-study feature sets
#   id_correlations_study4.csv — Study 4 feature x questionnaire r matrix
#   id_correlations_study5.csv — Study 5 feature x questionnaire r matrix
#   id_replication_study4_5.csv — replication comparison
#
# Assumes loaded: s1l, s1s, s2l, s2s, s3l, s3s, s4l, s4s,
#                 s5l, s5s, s5_long_breath, RESULTS_DIR
# ============================================================


source(file.path(ANALYSIS_DIR, "theme_bcat.R"))

# ── Helper: per-person OLS slope ──────────────────────────────

.pp_slope <- function(data, x_col, y_col, id_col = "id",
                       min_trials = 3) {
  data |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::group_map(function(g, key) {
      g <- tidyr::drop_na(g, dplyr::all_of(c(x_col, y_col)))
      if (nrow(g) < min_trials) return(NULL)
      x <- g[[x_col]]
      if (stats::sd(x, na.rm = TRUE) < 1e-10) return(NULL)
      tibble::tibble(
        id    = key[[id_col]],
        slope = stats::coef(stats::lm(
          stats::as.formula(paste(y_col, "~", x_col)),
          data = g))[[x_col]]
      )
    }, .keep = TRUE) |>
    purrr::compact() |>
    dplyr::bind_rows()
}

# ── Helper: per-person r (within-person correlation) ──────────

.pp_r <- function(data, x_col, y_col, id_col = "id",
                   min_trials = 5) {
  data |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::group_map(function(g, key) {
      g <- tidyr::drop_na(g, dplyr::all_of(c(x_col, y_col)))
      if (nrow(g) < min_trials) return(NULL)
      if (stats::sd(g[[x_col]], na.rm=TRUE) < 1e-10 ||
          stats::sd(g[[y_col]], na.rm=TRUE) < 1e-10) return(NULL)
      tibble::tibble(
        id = key[[id_col]],
        r  = stats::cor(g[[x_col]], g[[y_col]],
                        use = "complete.obs")
      )
    }, .keep = TRUE) |>
    purrr::compact() |>
    dplyr::bind_rows()
}

# ── Helper: extract features from long data ───────────────────

.extract_task_features <- function(long_data, summary_data,
                                    study_label,
                                    id_col      = "id",
                                    change_col  = "Change",
                                    arousal_col = "Arousal",
                                    conf_col    = "Confidence",
                                    acc_col     = "Accuracy") {

  # Coerce id to character in both datasets to avoid type mismatch
  long_data    <- dplyr::mutate(long_data,
                                 !!id_col := as.character(.data[[id_col]]))
  summary_data <- dplyr::mutate(summary_data,
                                 !!id_col := as.character(.data[[id_col]]))

  # Base: change trials only
  base <- long_data |>
    dplyr::filter(Direction %in% c("Faster","Slower"),
                  !is.na(.data[[acc_col]]),
                  !is.na(.data[[change_col]]),
                  !is.na(.data[[arousal_col]]))

  hits   <- dplyr::filter(base, .data[[acc_col]] == 1)
  misses <- dplyr::filter(base, .data[[acc_col]] == 0)
  acc_trials <- dplyr::filter(base, Direction == "Faster")
  dec_trials <- dplyr::filter(base, Direction == "Slower")

  message(sprintf("  [%s] N=%d participants", study_label,
                  dplyr::n_distinct(base[[id_col]])))

  # ── A. Task performance features ────────────────────────────

  # Miss rate
  miss_rate <- base |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::summarise(miss_rate = mean(.data[[acc_col]] == 0),
                     n_trials  = dplyr::n(), .groups = "drop")

  # Mean arousal, confidence by accuracy
  means_by_acc <- base |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::summarise(
      mean_arousal_overall = mean(.data[[arousal_col]], na.rm = TRUE),
      .groups = "drop")

  # ── B. Within-person slopes: Arousal ~ Change ────────────────

  # Overall hit/miss slopes
  slope_hit_overall  <- .pp_slope(hits,   change_col, arousal_col, id_col) |>
    dplyr::rename(slope_hit_overall  = slope)
  slope_miss_overall <- .pp_slope(misses, change_col, arousal_col, id_col) |>
    dplyr::rename(slope_miss_overall = slope)

  # Acceleration hit/miss slopes
  slope_hit_acc  <- .pp_slope(dplyr::filter(hits,   Direction=="Faster"),
                               change_col, arousal_col, id_col) |>
    dplyr::rename(slope_hit_acc = slope)
  slope_miss_acc <- .pp_slope(dplyr::filter(misses, Direction=="Faster"),
                               change_col, arousal_col, id_col) |>
    dplyr::rename(slope_miss_acc = slope)

  # Deceleration hit/miss slopes
  slope_hit_dec  <- .pp_slope(dplyr::filter(hits,   Direction=="Slower"),
                               change_col, arousal_col, id_col) |>
    dplyr::rename(slope_hit_dec = slope)
  slope_miss_dec <- .pp_slope(dplyr::filter(misses, Direction=="Slower"),
                               change_col, arousal_col, id_col) |>
    dplyr::rename(slope_miss_dec = slope)

  # ── C. Confidence-accuracy calibration ───────────────────────
  # Within-person r(confidence, accuracy): how well does
  # confidence track detection success?

  conf_acc_r <- if (!is.null(conf_col) && conf_col %in% names(base)) {
    .pp_r(base, conf_col, acc_col, id_col) |>
      dplyr::rename(conf_acc_r = r)
  } else {
    tibble::tibble(!!id_col := character(0), conf_acc_r = numeric(0))
  }

  # Within-person r(arousal, confidence)
  arousal_conf_r <- if (!is.null(conf_col) && conf_col %in% names(base)) {
    .pp_r(base, arousal_col, conf_col, id_col) |>
      dplyr::rename(arousal_conf_r = r)
  } else {
    tibble::tibble(!!id_col := character(0), arousal_conf_r = numeric(0))
  }

  # ── D. Arousal sensitivity: salience difference ───────────────
  # If Salience column available
  arousal_salience_diff <- if ("Salience" %in% names(base)) {
    base |>
      dplyr::group_by(.data[[id_col]], Salience) |>
      dplyr::summarise(m = mean(.data[[arousal_col]], na.rm=TRUE),
                       .groups="drop") |>
      tidyr::pivot_wider(names_from = Salience,
                         values_from = m,
                         names_prefix = "arousal_sal_") |>
      dplyr::mutate(arousal_salience_diff =
                      arousal_sal_High - arousal_sal_Low) |>
      dplyr::select(dplyr::all_of(id_col), arousal_salience_diff)
  } else {
    tibble::tibble(!!id_col := character(0),
                   arousal_salience_diff = numeric(0))
  }

  # ── E. Join all features ─────────────────────────────────────
  features <- miss_rate |>
    dplyr::left_join(means_by_acc,       by = id_col) |>
    dplyr::left_join(slope_hit_overall,  by = id_col) |>
    dplyr::left_join(slope_miss_overall, by = id_col) |>
    dplyr::left_join(slope_hit_acc,      by = id_col) |>
    dplyr::left_join(slope_miss_acc,     by = id_col) |>
    dplyr::left_join(slope_hit_dec,      by = id_col) |>
    dplyr::left_join(slope_miss_dec,     by = id_col) |>
    dplyr::left_join(conf_acc_r,         by = id_col) |>
    dplyr::left_join(arousal_conf_r,     by = id_col) |>
    dplyr::left_join(arousal_salience_diff, by = id_col)

  # Join summary-level features from summary_data
  summary_features <- summary_data |>
    dplyr::select(
      dplyr::all_of(id_col),
      # Thresholds
      dplyr::any_of(c("thresh_c1",
                       "thresh_Low_Faster","thresh_Low_Slower",
                       "thresh_High_Faster","thresh_High_Slower",
                       "thresh_ses1_Low_Faster","thresh_ses1_Low_Slower",
                       "thresh_ses1_High_Faster","thresh_ses1_High_Slower")),
      # Summary confidence and awareness
      dplyr::any_of(c("mean_Confidence","mean_Confidence_ses1",
                       "Awareness","Awareness_ses1",
                       "c_bias","c_bias_ses1",
                       "H_corr","H_corr_ses1",
                       "FA_corr","FA_corr_ses1",
                       "dprime","dprime_ses1")),
      # Questionnaires
      dplyr::any_of(c(
        # MAIA
        "MAIA_total","MAIA_Noticing","MAIA_NotDistracting",
        "MAIA_NotWorrying","MAIA_AttentionReg","MAIA_EmoAware",
        "MAIA_SelfReg","MAIA_BodyListen","MAIA_Trusting",
        # Body awareness
        "BARQ_total",
        # Stress/wellbeing
        "BIPS_total","PHQ4_Anxiety","PHQ4_Depression",
        "SPANE_Pos","SPANE_Neg","PSS_total","premood_stress",
        # Personality
        "MSES","MSES_selfdoubt","Alexithymia","BFI_Neuroticism",
        # Risk/reward
        "BARTtotalpoints","BIS","BAS_Drive","BAS_FunSeeking",
        # Affect
        "PANAS_Pos","PANAS_Neg"
      ))
    )

  # Compute mean threshold if condition-specific ones exist
  thresh_cols <- dplyr::intersect(
    c("thresh_Low_Faster","thresh_Low_Slower",
      "thresh_High_Faster","thresh_High_Slower",
      "thresh_ses1_Low_Faster","thresh_ses1_Low_Slower",
      "thresh_ses1_High_Faster","thresh_ses1_High_Slower"),
    names(summary_features))

  if (length(thresh_cols) > 0) {
    summary_features <- summary_features |>
      dplyr::mutate(
        mean_threshold = rowMeans(
          dplyr::across(dplyr::all_of(thresh_cols)), na.rm = TRUE),
        thresh_salience_diff = rowMeans(
          dplyr::across(dplyr::matches("High_")), na.rm = TRUE) -
          rowMeans(dplyr::across(dplyr::matches("Low_")), na.rm = TRUE),
        thresh_direction_diff = rowMeans(
          dplyr::across(dplyr::matches("_Slower")), na.rm = TRUE) -
          rowMeans(dplyr::across(dplyr::matches("_Faster")), na.rm = TRUE)
      )
  }

  features <- features |>
    dplyr::left_join(summary_features, by = id_col) |>
    dplyr::mutate(study = study_label)

  features
}


# ── Questionnaire columns for correlation analysis ────────────

qs_cols_s4 <- c(
  "MAIA_total","MAIA_Noticing","MAIA_NotDistracting","MAIA_NotWorrying",
  "MAIA_AttentionReg","MAIA_EmoAware","MAIA_SelfReg","MAIA_BodyListen",
  "MAIA_Trusting","BARQ_total","BIPS_total","PHQ4_Anxiety","PHQ4_Depression",
  "SPANE_Pos","SPANE_Neg","premood_stress")

qs_cols_s5 <- c(
  "MAIA_total","MAIA_Noticing","MAIA_NotDistracting","MAIA_NotWorrying",
  "MAIA_AttentionReg","MAIA_EmoAware","MAIA_SelfReg","MAIA_BodyListen",
  "MAIA_Trusting","BARQ_total","BIPS_total","PHQ4_Anxiety","PHQ4_Depression",
  "SPANE_Pos","SPANE_Neg","MSES_selfdoubt","Alexithymia")

task_feature_cols <- c(
  "miss_rate","mean_arousal_overall","mean_threshold",
  "thresh_salience_diff","thresh_direction_diff",
  "slope_hit_overall","slope_miss_overall",
  "slope_hit_acc","slope_miss_acc",
  "slope_hit_dec","slope_miss_dec",
  "conf_acc_r","arousal_conf_r","arousal_salience_diff",
  "mean_Confidence","Awareness","c_bias","H_corr","FA_corr","dprime")


# ── Extract features ──────────────────────────────────────────

message("\n========================================")
message("Individual Differences: Feature Extraction")
message("========================================")

# Study 4
message("\nStudy 4:")
feat_s4 <- .extract_task_features(
  s4l |> dplyr::filter(Group == "Breath"),
  s4s, "Study4")

# Study 5
message("\nStudy 5:")
feat_s5 <- .extract_task_features(
  s5_long_breath, s5s, "Study5",
  conf_col    = "Confidence",
  arousal_col = "Arousal")

# Study 1A (task features only — limited questionnaires)
message("\nStudy 1A:")
feat_s1 <- .extract_task_features(
  s1l |> dplyr::filter(Group == "TaskA"),
  s1s, "Study1A")

# Study 2
message("\nStudy 2:")
feat_s2 <- .extract_task_features(s2l, s2s, "Study2")

# Study 3
message("\nStudy 3:")
feat_s3 <- .extract_task_features(s3l, s3s, "Study3")

# Save all feature sets
for (feat_df in list(feat_s4, feat_s5, feat_s1, feat_s2, feat_s3)) {
  study_lbl <- feat_df$study[1]
  readr::write_csv(feat_df,
    file.path(RESULTS_DIR, paste0("id_features_",
           tolower(study_lbl), ".csv")))
  message(sprintf("Saved: id_features_%s.csv (%d participants, %d features)",
                  tolower(study_lbl), nrow(feat_df), ncol(feat_df)))
}


# ── Correlation matrices: task features x questionnaires ──────

.corr_matrix <- function(feat_df, task_cols, qs_cols,
                           study_label, min_n = 20) {

  # Retain only available columns
  task_avail <- dplyr::intersect(task_cols, names(feat_df))
  qs_avail   <- dplyr::intersect(qs_cols,   names(feat_df))

  results <- purrr::map_dfr(task_avail, function(tc) {
    purrr::map_dfr(qs_avail, function(qc) {
      x <- feat_df[[tc]]
      y <- feat_df[[qc]]
      complete <- !is.na(x) & !is.na(y)
      n <- sum(complete)
      if (n < min_n) return(NULL)
      ct <- tryCatch(
        stats::cor.test(x[complete], y[complete]),
        error = function(e) NULL)
      if (is.null(ct)) return(NULL)
      tibble::tibble(
        study        = study_label,
        task_feature = tc,
        questionnaire = qc,
        n            = n,
        r            = ct$estimate,
        t            = ct$statistic,
        p            = ct$p.value,
        ci_lower     = ct$conf.int[1],
        ci_upper     = ct$conf.int[2]
      )
    })
  })

  # BH correction within study
  if (nrow(results) > 0) {
    results <- results |>
      dplyr::mutate(p_bh = stats::p.adjust(p, method = "BH")) |>
      dplyr::arrange(p)
  }
  results
}

message("\n--- Correlation matrices ---")

corr_s4 <- .corr_matrix(feat_s4, task_feature_cols, qs_cols_s4, "Study4")
corr_s5 <- .corr_matrix(feat_s5, task_feature_cols, qs_cols_s5, "Study5")

readr::write_csv(corr_s4, file.path(RESULTS_DIR, "id_correlations_study4.csv"))
readr::write_csv(corr_s5, file.path(RESULTS_DIR, "id_correlations_study5.csv"))
message("Saved: id_correlations_study4.csv")
message("Saved: id_correlations_study5.csv")

# Print top correlations surviving BH correction
message("\nTop correlations surviving BH < .05 in Study 4:")
print(corr_s4 |>
  dplyr::filter(p_bh < .05) |>
  dplyr::mutate(across(c(r,p,p_bh), round, 3)) |>
  dplyr::select(task_feature, questionnaire, n, r, p, p_bh))

message("\nTop correlations surviving BH < .05 in Study 5:")
print(corr_s5 |>
  dplyr::filter(p_bh < .05) |>
  dplyr::mutate(across(c(r,p,p_bh), round, 3)) |>
  dplyr::select(task_feature, questionnaire, n, r, p, p_bh))


# ── Replication check: Study 4 → Study 5 ─────────────────────
# For each correlation found in Study 4, test whether it
# replicates (same direction, p < .05) in Study 5.

message("\n--- Study 4 → Study 5 replication check ---")

shared_qs <- dplyr::intersect(qs_cols_s4, qs_cols_s5)
shared_tf <- dplyr::intersect(task_feature_cols, task_feature_cols)

corr_s4_shared <- corr_s4 |>
  dplyr::filter(questionnaire %in% shared_qs)
corr_s5_shared <- corr_s5 |>
  dplyr::filter(questionnaire %in% shared_qs)

replication <- corr_s4_shared |>
  dplyr::inner_join(
    corr_s5_shared |>
      dplyr::select(task_feature, questionnaire,
                    r_s5 = r, p_s5 = p, p_bh_s5 = p_bh, n_s5 = n),
    by = c("task_feature","questionnaire")
  ) |>
  dplyr::rename(r_s4 = r, p_s4 = p, p_bh_s4 = p_bh, n_s4 = n) |>
  dplyr::mutate(
    same_direction  = sign(r_s4) == sign(r_s5),
    replicated      = same_direction & p_s5 < .05,
    sig_s4          = p_bh_s4 < .05,
    sig_s5          = p_bh_s5 < .05
  ) |>
  dplyr::arrange(p_s4)

readr::write_csv(replication,
                 file.path(RESULTS_DIR, "id_replication_study4_5.csv"))
message("Saved: id_replication_study4_5.csv")

# Summary
n_sig_s4    <- sum(replication$sig_s4, na.rm = TRUE)
n_replicated <- sum(replication$sig_s4 & replication$replicated,
                     na.rm = TRUE)
n_same_dir   <- sum(replication$sig_s4 & replication$same_direction,
                     na.rm = TRUE)

message(sprintf(
  "\nStudy 4 BH-significant correlations: %d", n_sig_s4))
message(sprintf(
  "Same direction in Study 5:            %d / %d (%.0f%%)",
  n_same_dir, n_sig_s4, 100*n_same_dir/max(n_sig_s4,1)))
message(sprintf(
  "Replicated (same dir + p<.05 in S5):  %d / %d (%.0f%%)",
  n_replicated, n_sig_s4, 100*n_replicated/max(n_sig_s4,1)))

message("\nTop replicated correlations:")
print(replication |>
  dplyr::filter(sig_s4, replicated) |>
  dplyr::mutate(across(c(r_s4, r_s5, p_s4, p_s5), round, 3)) |>
  dplyr::select(task_feature, questionnaire,
                r_s4, p_s4, r_s5, p_s5, n_s4, n_s5))
