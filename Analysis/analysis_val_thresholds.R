# ============================================================
# analysis_val_thresholds.R
# BCAT — Task Validation: Thresholds, d', ICC (Supplement S1.4)
#
# Corresponds to: Supplement S1.4 (H1 direction, H2 salience,
#                 H5 test-block d', H6 test-retest reliability)
#
# NOTE: analysis_val_detection.R must be sourced before this file.
#       table_validation.csv assembly uses s{N}_det objects created there.
#
# NOTE: s{N}_thresh, s{N}_conf, and s{N}_maia objects created here are
#       used downstream by analysis_maia.R. Run this file before
#       sourcing analysis_maia.R.
#
# Assumes in environment (set by MainAnalysis.R):
#   s1l, s1s, s2l, s2s, s3l, s3s, s4l, s4s, s4t, s5l, s5s, s5t
#   s1l_A_sorted, s5_long_breath
#   RESULTS_DIR
#   From analysis_val_detection.R: s{N}_det objects
#   Helpers from utils.R: make_threshold_long(), compute_retest_icc(),
#                         add_partial_r(), compute_test_dprime_3afc(),
#                         fmt_lmer(), fmt_cor()
#
# Creates in environment (used downstream by analysis_maia.R):
#   s1_thresh_A, s1_thresh_B, s4_thresh, s5_thresh
#   s1_conf_A, s1_conf_B, s2_conf, s3_conf, s4_conf, s5_conf
#   s1_maia_A, s1_maia_B, s2_maia, s3_maia, s4_maia, s5_maia
#   s4_thresh_long, s5_thresh_long
#   s4_test, s5_test
#
# Outputs:
#   table_reliability.csv
#   table_validation.csv
# ============================================================

message("\n========================================")
message("TASK VALIDATION: Thresholds, d', ICC (S1.4)")
message("========================================")


# ============================================================
# fit_confidence_models()
#
# Sensibility LMM: how strongly do confidence ratings track
# the objective change signal? Same hierarchical sequence as
# fit_detection_models() but on continuous Confidence.
#
# `data`: trial-level; canonical columns id, Confidence, Change.
# ============================================================
fit_confidence_models <- function(data, study_label = "") {

  data <- data |>
    dplyr::filter(!is.na(Confidence), !is.na(Change)) |>
    dplyr::mutate(Change2 = Change^2)

  cat(sprintf("\n[%s] Confidence models: N=%d trials, %d participants\n",
              study_label, nrow(data), dplyr::n_distinct(data$id)))

  m0  <- lmerTest::lmer(Confidence ~ 1               + (1 | id),
                    data = data, REML = FALSE)
  m1  <- lmerTest::lmer(Confidence ~ Change           + (1 | id),
                    data = data, REML = FALSE)
  m2  <- lmerTest::lmer(Confidence ~ Change + Change2 + (1 | id),
                    data = data, REML = FALSE)

  results <- list(
    null      = m0,
    linear    = m1,
    quadratic = m2,
    lrt_main  = anova(m0, m1, m2),
    study     = study_label
  )

  if ("Salience" %in% names(data)) {
    data <- data |>
      dplyr::mutate(Salience = factor(Salience, levels = c("Low", "High")))

    m3 <- lmerTest::lmer(
      Confidence ~ Change + Change2 + Salience + (1 | id),
      data = data, REML = FALSE
    )
    m4 <- lmerTest::lmer(
      Confidence ~ (Change + Change2) * Salience + (1 | id),
      data = data, REML = FALSE
    )
    results$salience_add <- m3
    results$salience_int <- m4
    results$lrt_salience <- anova(m2, m3, m4)
  }

  results$final <- update(m2, REML = TRUE)
  results
}


# ============================================================
# fit_maia_tests()
#
# MAIA dissociation tests (H3A, H3B).
#   H3A: Confidence ~ MAIA (expect positive r)
#   H3B: Threshold  ~ MAIA (expect null; Bayesian BF for H0)
#
# Optional partial correlation controlling for MSES_selfdoubt
# (Study 5 robustness: ruling out domain-general confidence bias).
#
# `data`:           one row per participant; must have MAIA_total_z,
#                   mean_Confidence, and at least one threshold column.
# `threshold_col`:  column name for H3B predictor.
# `confidence_col`: column name for H3A predictor (default "mean_Confidence").
# `control_var`:    optional column for partial correlation.
# ============================================================
fit_maia_tests <- function(data,
                           threshold_col  = "mean_Threshold",
                           confidence_col = "mean_Confidence",
                           control_var    = NULL,
                           study_label    = "") {

  data <- data |> dplyr::filter(!is.na(MAIA_total))

  if (!"MAIA_total_z" %in% names(data)) {
    warning(sprintf(
      "[%s] MAIA_total_z not found. Call standardise_maia() on summary data first. ",
      study_label
    ), "Computing within this call as fallback.")
    data <- standardise_maia(data, study_label = study_label)
  }

  results <- list(study = study_label)

  # ── H3A: Confidence ~ MAIA ──────────────────────────────────
  if (confidence_col %in% names(data)) {
    d3a <- data |> dplyr::filter(!is.na(.data[[confidence_col]]))

    results$H3A_freq <- cor.test(d3a[[confidence_col]], d3a$MAIA_total_z)
    results$H3A_n    <- nrow(d3a)

    cat(sprintf("\n[%s] H3A: %s ~ MAIA_z: %s\n",
                study_label, confidence_col,
                fmt_cor(results$H3A_freq)))

    if (!is.null(control_var) && control_var %in% names(d3a)) {
      d3a_pc <- d3a |>
        dplyr::select(conf = dplyr::all_of(confidence_col),
                      maia = MAIA_total_z,
                      ctrl = dplyr::all_of(control_var)) |>
        tidyr::drop_na()

      resid_conf <- residuals(lm(conf ~ ctrl, data = d3a_pc))
      resid_maia <- residuals(lm(maia ~ ctrl, data = d3a_pc))

      results$H3A_partial      <- cor.test(resid_conf, resid_maia)
      results$H3A_partial_n    <- nrow(d3a_pc)
      results$H3A_partial_ctrl <- control_var
      results$H3A_regression   <- lm(conf ~ maia + ctrl, data = d3a_pc)

      cat(sprintf("  Partial (controlling %s): %s  (N=%d)\n",
                  control_var,
                  fmt_cor(results$H3A_partial),
                  results$H3A_partial_n))
    }

    if ("Awareness" %in% names(d3a)) {
      d_aw <- d3a |> dplyr::filter(!is.na(Awareness))
      results$Awareness_MAIA   <- cor.test(d_aw$Awareness, d_aw$MAIA_total_z)
      results$Awareness_MAIA_n <- nrow(d_aw)
      cat(sprintf("  Awareness ~ MAIA_z: %s\n",
                  fmt_cor(results$Awareness_MAIA)))

      if (!is.null(control_var) && control_var %in% names(d_aw)) {
        d_aw_pc <- d_aw |>
          dplyr::select(Awareness, maia = MAIA_total_z,
                        ctrl = dplyr::all_of(control_var)) |>
          tidyr::drop_na()
        resid_aw   <- residuals(lm(Awareness ~ ctrl, data = d_aw_pc))
        resid_maia <- residuals(lm(maia       ~ ctrl, data = d_aw_pc))
        results$Awareness_MAIA_partial   <- cor.test(resid_aw, resid_maia)
        results$Awareness_MAIA_partial_n <- nrow(d_aw_pc)
      }
    }
  }

  # ── H3B: Threshold ~ MAIA ───────────────────────────────────
  if (threshold_col %in% names(data)) {
    d3b <- data |> dplyr::filter(!is.na(.data[[threshold_col]]))

    results$H3B_freq <- cor.test(d3b[[threshold_col]], d3b$MAIA_total_z)
    results$H3B_n    <- nrow(d3b)

    bf_result <- tryCatch(
      BayesFactor::correlationBF(
        d3b[[threshold_col]],
        d3b$MAIA_total,
        iterations = 10000
      ),
      error = function(e) NULL
    )
    if (!is.null(bf_result)) {
      results$H3B_BF   <- BayesFactor::extractBF(bf_result)$bf
      results$H3B_BF01 <- 1 / results$H3B_BF

      cat(sprintf("[%s] H3B: %s ~ MAIA_z: %s  BF01=%.2f\n",
                  study_label, threshold_col,
                  fmt_cor(results$H3B_freq),
                  results$H3B_BF01))
    }
  }

  results
}


# ============================================================
# fit_threshold_models()
#
# Threshold LMMs for H1 (Direction) and H2 (Salience).
# Operates on long-format threshold data: one row per
# participant per Salience x Direction condition (Studies 4, 5)
# or per task (Study 1: one row per TaskA/TaskB per participant).
#
# `has_salience`: fit H2 and interaction (FALSE for Studies 1, 2).
# `has_group`:    include Group as additional fixed effect.
# ============================================================
fit_threshold_models <- function(data,
                                 has_salience = TRUE,
                                 has_group    = FALSE,
                                 study_label  = "") {

  data <- data |>
    dplyr::filter(!is.na(Threshold)) |>
    dplyr::mutate(
      Direction = factor(Direction, levels = c("Faster", "Slower"))
    )

  cat(sprintf("\n[%s] Threshold models: N=%d observations, %d participants\n",
              study_label, nrow(data), dplyr::n_distinct(data$id)))

  results <- list(study = study_label)

  m_dir <- lmerTest::lmer(Threshold ~ Direction + (1 | id),
                          data = data, REML = TRUE)
  results$H1     <- m_dir
  results$H1_emm <- emmeans::emmeans(m_dir, ~ Direction)

  cat(sprintf("  H1 Direction: %s\n", fmt_lmer(m_dir, "DirectionSlower")))

  if (has_salience) {
    data <- data |>
      dplyr::mutate(Salience = factor(Salience, levels = c("Low", "High")))

    m_sal <- lmerTest::lmer(Threshold ~ Salience + (1 | id),
                            data = data, REML = TRUE)
    results$H2     <- m_sal
    results$H2_emm <- emmeans::emmeans(m_sal, ~ Salience)

    cat(sprintf("  H2 Salience: %s\n", fmt_lmer(m_sal, "SalienceHigh")))

    m_add <- lmerTest::lmer(Threshold ~ Direction + Salience       + (1 | id),
                            data = data, REML = TRUE)
    m_int <- lmerTest::lmer(Threshold ~ Direction * Salience       + (1 | id),
                            data = data, REML = TRUE)
    results$additive    <- m_add
    results$interaction <- m_int
    results$lrt_dir_sal <- anova(
      update(m_dir, REML = FALSE),
      update(m_sal, REML = FALSE),
      update(m_add, REML = FALSE),
      update(m_int, REML = FALSE)
    )

    int_cf  <- summary(m_int)$coefficients
    int_row <- grep("Direction.*Salience|Salience.*Direction",
                    rownames(int_cf), value = TRUE)
    if (length(int_row) > 0)
      results$interaction_p <- int_cf[int_row[1], "Pr(>|t|)"]
  }

  if (has_group && "Group" %in% names(data)) {
    m_grp <- lmerTest::lmer(
      Threshold ~ Direction + Salience + Group + (1 | id),
      data = data, REML = TRUE
    )
    results$with_group <- m_grp
  }

  results
}


# ============================================================
# fit_test_accuracy()
#
# Test-phase sensitivity check (H5). Predicts binary Accuracy
# from Salience and Direction on fixed-level test trials.
# Uses glmer (binomial). Filter to Condition == "breath" upstream.
# ============================================================
fit_test_accuracy <- function(data,
                               include_group = FALSE,
                               study_label   = "") {

  data <- data |>
    dplyr::filter(!is.na(Accuracy)) |>
    dplyr::mutate(
      Accuracy  = as.integer(Accuracy),
      Salience  = factor(Salience,  levels = c("Low",    "High")),
      Direction = factor(Direction, levels = c("Faster", "Slower"))
    )

  cat(sprintf("\n[%s] Test-phase accuracy: N=%d trials, %d participants\n",
              study_label, nrow(data), dplyr::n_distinct(data$id)))

  glmer_ctrl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  m0  <- lme4::glmer(Accuracy ~ 1                    + (1 | id),
                      data = data, family = binomial, control = glmer_ctrl)
  m_s <- lme4::glmer(Accuracy ~ Salience             + (1 | id),
                      data = data, family = binomial, control = glmer_ctrl)
  m_d <- lme4::glmer(Accuracy ~ Salience + Direction + (1 | id),
                      data = data, family = binomial, control = glmer_ctrl)
  m_i <- lme4::glmer(Accuracy ~ Salience * Direction + (1 | id),
                      data = data, family = binomial, control = glmer_ctrl)

  results <- list(
    H5_null = m0,
    H5_sal  = m_s,
    H5_add  = m_d,
    H5_int  = m_i,
    lrt     = anova(m0, m_s, m_d, m_i),
    study   = study_label
  )

  cat(sprintf("  H5 Salience effect: OR=%.2f, p=%.4f\n",
              exp(lme4::fixef(m_s)["SalienceHigh"]),
              summary(m_s)$coefficients["SalienceHigh", "Pr(>|z|)"]))

  if (include_group && "Group" %in% names(data)) {
    m_grp <- lme4::glmer(
      Accuracy ~ Salience + Direction + Group + (1 | id),
      data = data, family = binomial, control = glmer_ctrl
    )
    results$H5_group <- m_grp
  }

  results
}


# ── Run models ────────────────────────────────────────────────

# Confidence models (all studies)
s1_conf_A <- fit_confidence_models(
  dplyr::filter(s1l, Group == "TaskA"), study_label = "Study1_TaskA")
s1_conf_B <- fit_confidence_models(
  dplyr::filter(s1l, Group == "TaskB"), study_label = "Study1_TaskB")
s2_conf   <- fit_confidence_models(s2l, study_label = "Study2")
s3_conf   <- fit_confidence_models(s3l, study_label = "Study3")
s4_conf   <- fit_confidence_models(s4l, study_label = "Study4")
s5_conf   <- fit_confidence_models(s5l, study_label = "Study5")

# MAIA tests (all studies)
s1_maia_A <- fit_maia_tests(
  s1s |> dplyr::rename(mean_Confidence = mean_Confidence_TaskA) |>
    dplyr::mutate(mean_Threshold = thresh_TaskA, Awareness = Awareness_TaskA),
  threshold_col = "mean_Threshold", study_label = "Study1_TaskA")
s1_maia_B <- fit_maia_tests(
  s1s |> dplyr::rename(mean_Confidence = mean_Confidence_TaskB) |>
    dplyr::mutate(mean_Threshold = thresh_TaskB, Awareness = Awareness_TaskB),
  threshold_col = "mean_Threshold", study_label = "Study1_TaskB")
s2_maia <- fit_maia_tests(
  s2s |> dplyr::mutate(mean_Threshold = thresh_c1),
  threshold_col = "mean_Threshold", study_label = "Study2")
s3_maia <- fit_maia_tests(
  s3s |> dplyr::rename(mean_Confidence = mean_Confidence_run1,
                        mean_Threshold  = mean_Accuracy),
  threshold_col = "mean_Threshold", study_label = "Study3")

# Threshold models: Studies 4 and 5 (H1, H2)
s4_thresh_long <- make_threshold_long(s4s, study = 4)
s4_thresh <- fit_threshold_models(
  s4_thresh_long, has_salience = TRUE, has_group = TRUE, study_label = "Study4")
s4_maia <- fit_maia_tests(
  s4s |> dplyr::mutate(
    mean_Threshold  = rowMeans(dplyr::across(dplyr::matches("^thresh_")), na.rm = TRUE),
    mean_Confidence = mean_Confidence),
  threshold_col = "mean_Threshold", study_label = "Study4")

s5s <- s5s |>
  dplyr::mutate(
    overall_mean_threshold = rowMeans(
      dplyr::across(dplyr::matches("^thresh_ses1_")), na.rm = TRUE)
  )
s5_thresh_long <- make_threshold_long(s5s, study = 5)
s5_thresh <- fit_threshold_models(
  dplyr::filter(s5_thresh_long, Direction %in% c("Faster", "Slower")),
  has_salience = TRUE, has_group = TRUE, study_label = "Study5")
s5_maia <- fit_maia_tests(
  s5s,
  threshold_col  = "overall_mean_threshold",
  confidence_col = "mean_Confidence",
  control_var    = "MSES_selfdoubt",
  study_label    = "Study5")

# Study 1 threshold direction models
s1_thresh_A <- fit_threshold_models(
  s1l |> dplyr::filter(Group == "TaskA") |>
    dplyr::rename(Threshold = TotalRateChange) |>
    dplyr::filter(Direction %in% c("Faster", "Slower")),
  has_salience = FALSE, study_label = "Study1_TaskA")
s1_thresh_B <- fit_threshold_models(
  s1l |> dplyr::filter(Group == "TaskB", Accuracy == 1,
                        Direction != "NoChange", DetectedEarly == TRUE) |>
    dplyr::mutate(Threshold = abs(as.numeric(DetectedChange) - 1)) |>
    dplyr::filter(Direction %in% c("Faster", "Slower")),
  has_salience = FALSE, study_label = "Study1_TaskB")


# ── Test-block accuracy GLMM (H5; Studies 4 and 5) ───────────
# NOTE: test_block_accuracy.R now generates the full direction- and
# session-stratified d' table (test_block_dprime_by_group.csv) used
# for Table ST5. The pooled d' by salience (table_test_dprime.csv)
# is no longer needed and has been removed. s4_test / s5_test are
# kept below because their H5_sal model objects feed table_validation
# and analysis_maia.R.

s4_test <- fit_test_accuracy(
  dplyr::filter(s4t, !is.na(Salience)),
  include_group = TRUE, study_label = "Study4")
s5_test <- fit_test_accuracy(
  dplyr::filter(s5t, Condition == "breath", !is.na(Salience)),
  include_group = TRUE, study_label = "Study5")


# ── Test-retest ICC (H6; Study 5 only) ───────────────────────
#
# Pre-registered criterion: ICC >= .70. Not met at condition level
# (range .24-.59). Aggregate r (breath-first r = .673,
# visual-first r = .318). See pre-registration deviations note.

thresh_conditions <- c("High_Faster", "High_Slower", "Low_Faster", "Low_Slower")

s5s_breath <- dplyr::filter(s5s, Group == "breath")
s5s_visual <- dplyr::filter(s5s, Group == "visual")

s5_icc_breath <- purrr::map_dfr(thresh_conditions, function(cond) {
  res <- compute_retest_icc(
    s5s_breath,
    ses1_col = paste0("thresh_ses1_", cond),
    ses2_col = paste0("thresh_ses2_", cond),
    label    = cond
  )
  tibble::tibble(
    condition = cond, group = "Breath-first (breath-breath)",
    n = res$n, icc = res$icc$value,
    icc_lower = res$icc$lbound, icc_upper = res$icc$ubound
  )
})

s5_icc_visual <- purrr::map_dfr(thresh_conditions, function(cond) {
  res <- compute_retest_icc(
    s5s_visual,
    ses1_col = paste0("thresh_ses1_", cond),
    ses2_col = paste0("thresh_ses2_", cond),
    label    = cond
  )
  tibble::tibble(
    condition = cond, group = "Visual-first (visual-breath, cross-condition)",
    n = res$n, icc = res$icc$value,
    icc_lower = res$icc$lbound, icc_upper = res$icc$ubound
  )
})

s5_agg_breath <- s5s_breath |>
  dplyr::mutate(
    avg_ses1 = rowMeans(dplyr::across(dplyr::matches("thresh_ses1")), na.rm = TRUE),
    avg_ses2 = rowMeans(dplyr::across(dplyr::matches("thresh_ses2")), na.rm = TRUE)
  ) |> dplyr::filter(!is.na(avg_ses1), !is.na(avg_ses2))
s5_agg_visual <- s5s_visual |>
  dplyr::mutate(
    avg_ses1 = rowMeans(dplyr::across(dplyr::matches("thresh_ses1")), na.rm = TRUE),
    avg_ses2 = rowMeans(dplyr::across(dplyr::matches("thresh_ses2")), na.rm = TRUE)
  ) |> dplyr::filter(!is.na(avg_ses1), !is.na(avg_ses2))

cat(sprintf(
  "\nH6 aggregate r (breath-first): r = %.3f (N=%d)\n",
  cor(s5_agg_breath$avg_ses1, s5_agg_breath$avg_ses2),
  nrow(s5_agg_breath)))
cat(sprintf(
  "H6 aggregate r (visual-first): r = %.3f (N=%d)\n",
  cor(s5_agg_visual$avg_ses1, s5_agg_visual$avg_ses2),
  nrow(s5_agg_visual)))

table_reliability <- dplyr::bind_rows(s5_icc_breath, s5_icc_visual) |>
  dplyr::mutate(study = "Study5")
readr::write_csv(table_reliability,
                 file.path(RESULTS_DIR, "table_reliability.csv"))
message("Saved: table_reliability.csv")


# ── Assemble cross-study validation table ────────────────────
#
# NOTE: uses s{N}_det objects from analysis_val_detection.R.
#       That file must be sourced first (guaranteed by MainAnalysis.R order).

.extract_coef <- function(model, term_pattern, p_col = "Pr(>|t|)") {
  if (is.null(model)) return(c(b = NA_real_, se = NA_real_, p = NA_real_))
  cf <- tryCatch(summary(model)$coefficients, error = function(e) NULL)
  if (is.null(cf)) return(c(b = NA_real_, se = NA_real_, p = NA_real_))
  row_nm <- grep(term_pattern, rownames(cf), value = TRUE)[1]
  if (is.na(row_nm)) return(c(b = NA_real_, se = NA_real_, p = NA_real_))
  c(b  = cf[row_nm, "Estimate"],
    se = cf[row_nm, "Std. Error"],
    p  = cf[row_nm, p_col])
}

table_validation <- tibble::tibble(
  study = c("Study1A", "Study1B", "Study2", "Study4", "Study5"),
  H1_b  = c(.extract_coef(s1_thresh_A$H1, "DirectionSlower")["b"],
             .extract_coef(s1_thresh_B$H1, "DirectionSlower")["b"],
             NA_real_,
             .extract_coef(s4_thresh$H1,   "DirectionSlower")["b"],
             .extract_coef(s5_thresh$H1,   "DirectionSlower")["b"]),
  H1_SE = c(.extract_coef(s1_thresh_A$H1, "DirectionSlower")["se"],
             .extract_coef(s1_thresh_B$H1, "DirectionSlower")["se"],
             NA_real_,
             .extract_coef(s4_thresh$H1,   "DirectionSlower")["se"],
             .extract_coef(s5_thresh$H1,   "DirectionSlower")["se"]),
  H1_p  = c(.extract_coef(s1_thresh_A$H1, "DirectionSlower")["p"],
             .extract_coef(s1_thresh_B$H1, "DirectionSlower")["p"],
             NA_real_,
             .extract_coef(s4_thresh$H1,   "DirectionSlower")["p"],
             .extract_coef(s5_thresh$H1,   "DirectionSlower")["p"]),
  H2_b  = c(NA_real_, NA_real_, NA_real_,
             .extract_coef(s4_thresh$H2, "SalienceHigh")["b"],
             .extract_coef(s5_thresh$H2, "SalienceHigh")["b"]),
  H2_SE = c(NA_real_, NA_real_, NA_real_,
             .extract_coef(s4_thresh$H2, "SalienceHigh")["se"],
             .extract_coef(s5_thresh$H2, "SalienceHigh")["se"]),
  H2_p  = c(NA_real_, NA_real_, NA_real_,
             .extract_coef(s4_thresh$H2, "SalienceHigh")["p"],
             .extract_coef(s5_thresh$H2, "SalienceHigh")["p"]),
  H5_b  = c(NA_real_, NA_real_, NA_real_,
             .extract_coef(s4_test$H5_sal, "SalienceHigh", p_col = "Pr(>|z|)")["b"],
             .extract_coef(s5_test$H5_sal, "SalienceHigh", p_col = "Pr(>|z|)")["b"]),
  H5_SE = c(NA_real_, NA_real_, NA_real_,
             .extract_coef(s4_test$H5_sal, "SalienceHigh", p_col = "Pr(>|z|)")["se"],
             .extract_coef(s5_test$H5_sal, "SalienceHigh", p_col = "Pr(>|z|)")["se"]),
  H5_p  = c(NA_real_, NA_real_, NA_real_,
             .extract_coef(s4_test$H5_sal, "SalienceHigh", p_col = "Pr(>|z|)")["p"],
             .extract_coef(s5_test$H5_sal, "SalienceHigh", p_col = "Pr(>|z|)")["p"])
)
table_validation <- add_partial_r(table_validation, c("H1", "H2", "H5"))

readr::write_csv(table_validation,
                 file.path(RESULTS_DIR, "table_validation.csv"))
message("Saved: table_validation.csv")
