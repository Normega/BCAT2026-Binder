# ============================================================
# analysis_s7_maia_selfesteem.R
# Study 5: MAIA Specificity — Controlling for Self-Esteem
#
# Tests whether MAIA predicts interoceptive confidence and
# detection threshold over and above general self-esteem /
# self-doubt (MSES total and MSES_selfdoubt subscale).
#
# Argument: if MAIA-confidence correlations reflect genuine
# interoceptive sensibility rather than general self-confidence
# or trait self-doubt, they should survive controlling for MSES.
#
# Sections:
#   s7.1  Data preparation — z-scores, session-long pivot
#   s7.2  Partial correlations: MAIA → confidence (per session)
#   s7.3  Partial correlations: MAIA → threshold (overall + 4 conditions)
#   s7.4  Multilevel model: confidence ~ MAIA + MSES + condition
#           (a) visual-first group only — clean within-person design
#           (b) full sample — session + condition as covariates
#   s7.5  Multilevel model: threshold ~ MAIA + MSES + Salience + Direction
#
# DESIGN NOTE:
#   In Study 5, Session 2 is always breath (for both groups).
#   Session 1 is breath for breath-first, visual for visual-first.
#   Condition and Session are therefore partially confounded in the
#   full sample. The visual-first subgroup (s7.4a) provides the
#   cleanest within-person test of a condition × MAIA interaction
#   because it has both conditions across sessions.
#
# Sources after analysis_maia.R (requires s5s with MAIA_total_z,
# s5_thresh_long from analysis_val_thresholds.R).
# Adds MSES_z and MSES_selfdoubt_z to s5s locally.
#
# Outputs:
#   s7_partial_correlations.csv
#   s7_threshold_partial_correlations.csv
#   s7_multilevel_confidence_results.csv
#   s7_multilevel_threshold_results.csv
# ============================================================

message("\n========================================")
message("s7: MAIA SPECIFICITY — SELF-ESTEEM CONTROLS")
message("========================================")


# ============================================================
# s7.1  Data Preparation
# ============================================================

message("\n--- s7.1: Data preparation ---")

# Create s5_thresh_long if not already in environment
# (normally created by analysis_val_thresholds.R)
if (!exists("s5_thresh_long")) {
  message("  s5_thresh_long not found -- creating from s5s")
  s5_thresh_long <- make_threshold_long(s5s, study = 5)
}

ctrl <- lme4::lmerControl(optimizer = "bobyqa")

# Add MSES z-scores to s5s (MAIA_total_z already created by
# standardise_maia() in analysis_val_thresholds.R)
s5_se <- s5s |>
  dplyr::mutate(
    MSES_total_z     = as.numeric(scale(MSES)),
    MSES_selfdoubt_z = as.numeric(scale(MSES_selfdoubt)),
    # Mean threshold across all 4 staircase conditions per session
    # Using breath session only: ses1 for breath-first, ses2 for visual-first
    breath_ses = dplyr::if_else(Condition_ses1 == "breath", "ses1", "ses2"),
    thresh_breath_Low_Faster  = dplyr::if_else(
      breath_ses == "ses1", thresh_ses1_Low_Faster,  thresh_ses2_Low_Faster),
    thresh_breath_Low_Slower  = dplyr::if_else(
      breath_ses == "ses1", thresh_ses1_Low_Slower,  thresh_ses2_Low_Slower),
    thresh_breath_High_Faster = dplyr::if_else(
      breath_ses == "ses1", thresh_ses1_High_Faster, thresh_ses2_High_Faster),
    thresh_breath_High_Slower = dplyr::if_else(
      breath_ses == "ses1", thresh_ses1_High_Slower, thresh_ses2_High_Slower),
    thresh_breath_mean = rowMeans(
      dplyr::pick(thresh_breath_Low_Faster, thresh_breath_Low_Slower,
                  thresh_breath_High_Faster, thresh_breath_High_Slower),
      na.rm = TRUE),
    # Confidence per session
    conf_ses1 = mean_Confidence_ses1,
    conf_ses2 = mean_Confidence_ses2
  )

cat(sprintf("Study 5 N: %d total, %d with MSES, %d with MAIA_total_z\n",
            nrow(s5_se),
            sum(!is.na(s5_se$MSES)),
            sum(!is.na(s5_se$MAIA_total_z))))


# ============================================================
# s7.2  Partial Correlations: MAIA → Confidence
#
# For each session, test:
#   (a) zero-order: r(MAIA, Confidence)
#   (b) partial controlling for MSES_total
#   (c) partial controlling for MSES_selfdoubt
#   (d) partial controlling for both MSES measures
#
# Separately for MAIA_total and MAIA subscales most relevant
# to sensibility (BodyListen, Noticing) since these are
# the theoretically motivated predictors.
# ============================================================

message("\n--- s7.2: Partial correlations: MAIA → Confidence ---")

# Helper: run zero-order and partial correlations for one predictor
.run_partial_conf <- function(data, predictor, outcome, controls, label) {
  dat <- data |>
    dplyr::select(dplyr::all_of(c(predictor, outcome, unlist(controls)))) |>
    tidyr::drop_na()
  n <- nrow(dat)
  if (n < 10) {
    message(sprintf("  [SKIP] %s: too few complete cases (N=%d)", label, n))
    return(NULL)
  }

  # Zero-order
  r0 <- cor.test(dat[[predictor]], dat[[outcome]])

  # Partial correlations
  results <- list(tibble::tibble(
    label       = label,
    predictor   = predictor,
    outcome     = outcome,
    control     = "none",
    r           = r0$estimate,
    t_stat      = r0$statistic,
    df          = r0$parameter,
    p           = r0$p.value,
    ci_lo       = r0$conf.int[1],
    ci_hi       = r0$conf.int[2],
    n           = n
  ))

  for (ctrl_set in controls) {
    z_mat <- as.matrix(dat[, ctrl_set, drop = FALSE])
    pr    <- ppcor::pcor.test(dat[[predictor]], dat[[outcome]], z_mat)
    results <- c(results, list(tibble::tibble(
      label       = label,
      predictor   = predictor,
      outcome     = outcome,
      control     = ctrl_set,
      r           = pr$estimate,
      t_stat      = pr$statistic,
      df          = pr$n - 2 - pr$gp,
      p           = pr$p.value,
      ci_lo       = NA_real_,
      ci_hi       = NA_real_,
      n           = pr$n
    )))
  }

  dplyr::bind_rows(results)
}

# Control sets to test
control_sets <- list(
  "MSES_total_z",
  "MSES_selfdoubt_z",
  c("MSES_total_z", "MSES_selfdoubt_z")
)

# Predictors of interest
maia_preds <- c("MAIA_total_z", "MAIA_BodyListen_z", "MAIA_Noticing_z")

# Run across sessions and predictors
conf_partial_rows <- list()
for (ses_label in c("ses1", "ses2")) {
  conf_col   <- paste0("conf_", ses_label)
  ses_data   <- s5_se |> dplyr::filter(!is.na(.data[[conf_col]]))
  
  ses_data <- ses_data |>
    dplyr::mutate(
      session_condition = if (ses_label == "ses1") Condition_ses1 else "breath"
    )

  for (pred in maia_preds) {
    if (!pred %in% names(ses_data)) next
    row <- .run_partial_conf(
      data      = ses_data,
      predictor = pred,
      outcome   = conf_col,
      controls  = control_sets,
      label     = ses_label
    )
    conf_partial_rows <- c(conf_partial_rows, list(row))
  }
}

s7_conf_partials <- dplyr::bind_rows(conf_partial_rows)
print(s7_conf_partials |>
        dplyr::filter(predictor == "MAIA_total_z") |>
        dplyr::select(label, control, r, p, n))


# ============================================================
# s7.3  Partial Correlations: MAIA → Threshold
#
# Outcome: breath session thresholds (overall mean + 4 conditions)
# ============================================================

message("\n--- s7.3: Partial correlations: MAIA → Threshold ---")

thresh_outcomes <- c(
  "thresh_breath_mean",
  "thresh_breath_Low_Faster",
  "thresh_breath_Low_Slower",
  "thresh_breath_High_Faster",
  "thresh_breath_High_Slower"
)

thresh_partial_rows <- list()
for (thresh_col in thresh_outcomes) {
  for (pred in maia_preds) {
    if (!pred %in% names(s5_se)) next
    row <- .run_partial_conf(
      data      = s5_se,
      predictor = pred,
      outcome   = thresh_col,
      controls  = control_sets,
      label     = thresh_col
    )
    thresh_partial_rows <- c(thresh_partial_rows, list(row))
  }
}

s7_thresh_partials <- dplyr::bind_rows(thresh_partial_rows)
cat("\nMAIA_total_z → thresh_breath_mean (all control sets):\n")
print(s7_thresh_partials |>
        dplyr::filter(predictor == "MAIA_total_z",
                      outcome   == "thresh_breath_mean") |>
        dplyr::select(control, r, p, n))


# ============================================================
# s7.4  Multilevel Model: Confidence ~ MAIA + MSES + Condition
#
# Pivots to long: one row per participant × session.
# Confidence is the session-level mean.
#
# (a) Visual-first group only
#     ses1 = visual, ses2 = breath → clean within-person contrast.
#     Key term: MAIA_total_z:Condition
#     Interpretation: does MAIA predict confidence more strongly
#     during breath-focused than visual-tracking sessions?
#
# (b) Full sample
#     Session and condition partially confounded (ses2 always breath).
#     Include both as covariates; interpret Condition effect cautiously.
# ============================================================

message("\n--- s7.4: Multilevel models: Confidence ~ MAIA + MSES + Condition ---")

# Build session-long data
conf_long <- s5_se |>
  dplyr::select(id, Condition_ses1, completed_ses2,
                conf_ses1, conf_ses2,
                MAIA_total_z, MAIA_BodyListen_z, MAIA_Noticing_z,
                MSES_total_z, MSES_selfdoubt_z) |>
  tidyr::pivot_longer(
    cols      = c(conf_ses1, conf_ses2),
    names_to  = "ses",
    names_prefix = "conf_",
    values_to = "Confidence"
  ) |>
  dplyr::filter(!is.na(Confidence)) |>
  dplyr::mutate(
    # Condition for each session
    Condition = dplyr::case_when(
      ses == "ses1" ~ Condition_ses1,
      ses == "ses2" ~ "breath"
    ),
    Condition = factor(Condition, levels = c("visual", "breath")),
    ses       = factor(ses, levels = c("ses1", "ses2"))
  )

cat(sprintf("\nConf long data: %d rows, %d participants\n",
            nrow(conf_long), dplyr::n_distinct(conf_long$id)))
cat("Condition × session cross-tab:\n")
print(table(conf_long$ses, conf_long$Condition))

# ── s7.4a: Visual-first group (within-person, clean design) ──────────

message("\n  s7.4a: Visual-first group only")

conf_vf <- conf_long |>
  dplyr::filter(Condition_ses1 == "visual",
                completed_ses2 == TRUE)

cat(sprintf("  Visual-first N: %d participants, %d observations\n",
            dplyr::n_distinct(conf_vf$id), nrow(conf_vf)))

# Complete cases for all predictors used in any model in this section
conf_vf_cc <- conf_vf |>
  dplyr::filter(
    !is.na(Confidence), !is.na(MAIA_total_z),
    !is.na(MSES_total_z), !is.na(MSES_selfdoubt_z)
  )

cat(sprintf("  Visual-first N complete: %d participants, %d observations\n",
            dplyr::n_distinct(conf_vf_cc$id), nrow(conf_vf_cc)))

# Base model
m_conf_vf_base <- lmerTest::lmer(
  Confidence ~ Condition + ses + (1 | id),
  data = conf_vf_cc, REML = FALSE, control = ctrl
)

# MAIA main effect
m_conf_vf_maia <- lmerTest::lmer(
  Confidence ~ MAIA_total_z + Condition + ses + (1 | id),
  data = conf_vf_cc, REML = FALSE, control = ctrl
)

# MSES controls
m_conf_vf_mses <- lmerTest::lmer(
  Confidence ~ MAIA_total_z + MSES_total_z + MSES_selfdoubt_z +
    Condition + ses + (1 | id),
  data = conf_vf_cc, REML = FALSE, control = ctrl
)

# MAIA × Condition interaction (key test)
m_conf_vf_int <- lmerTest::lmer(
  Confidence ~ MAIA_total_z * Condition + MSES_total_z + MSES_selfdoubt_z +
    ses + (1 | id),
  data = conf_vf_cc, REML = FALSE, control = ctrl
)

cat("\n  Fixed effects (visual-first, MAIA × Condition model):\n")
print(round(summary(m_conf_vf_int)$coefficients, 4))

lrt_vf <- anova(m_conf_vf_base, m_conf_vf_maia, m_conf_vf_mses, m_conf_vf_int)
cat("\n  LRT sequence (visual-first):\n")
print(lrt_vf)

# ── s7.4b: Full sample ───────────────────────────────────────────────

message("\n  s7.4b: Full sample")

conf_full <- conf_long |>
  dplyr::filter(!is.na(MAIA_total_z),
                !is.na(MSES_total_z))

cat(sprintf("  Full sample N: %d participants, %d observations\n",
            dplyr::n_distinct(conf_full$id), nrow(conf_full)))

conf_full_cc <- conf_full |>
  dplyr::filter(
    !is.na(Confidence), !is.na(MAIA_total_z),
    !is.na(MSES_total_z), !is.na(MSES_selfdoubt_z)
  )

cat(sprintf("  Full sample N complete: %d participants, %d observations\n",
            dplyr::n_distinct(conf_full_cc$id), nrow(conf_full_cc)))


m_conf_full_base <- lmerTest::lmer(
  Confidence ~ Condition + ses + (1 | id),
  data = conf_full_cc, REML = FALSE, control = ctrl
)

m_conf_full_maia <- lmerTest::lmer(
  Confidence ~ MAIA_total_z + Condition + ses + (1 | id),
  data = conf_full_cc, REML = FALSE, control = ctrl
)

m_conf_full_mses <- lmerTest::lmer(
  Confidence ~ MAIA_total_z + MSES_total_z + MSES_selfdoubt_z +
    Condition + ses + (1 | id),
  data = conf_full_cc, REML = FALSE, control = ctrl
)

m_conf_full_int <- lmerTest::lmer(
  Confidence ~ MAIA_total_z * Condition + MSES_total_z + MSES_selfdoubt_z +
    ses + (1 | id),
  data = conf_full_cc, REML = FALSE, control = ctrl
)

cat("\n  Fixed effects (full sample, MAIA × Condition model):\n")
print(round(summary(m_conf_full_int)$coefficients, 4))

lrt_full <- anova(m_conf_full_base, m_conf_full_maia,
                  m_conf_full_mses, m_conf_full_int)
cat("\n  LRT sequence (full sample):\n")
print(lrt_full)

# Extract key coefficients for output table
.extract_mlm <- function(model, model_label, group_label) {
  cf <- summary(model)$coefficients
  broom.mixed::tidy(model, effects = "fixed") |>
    dplyr::mutate(
      model_label = model_label,
      group_label = group_label,
      partial_r   = statistic / sqrt(statistic^2 + df)
    )
}

s7_conf_mlm <- dplyr::bind_rows(
  .extract_mlm(m_conf_vf_int,   "MAIA*Condition + MSES", "visual-first"),
  .extract_mlm(m_conf_full_int, "MAIA*Condition + MSES", "full-sample")
)


# ============================================================
# s7.5  Multilevel Model: Threshold ~ MAIA + MSES + Salience + Direction
#
# Uses s5_thresh_long (from analysis_val_thresholds.R), restricted
# to breath sessions. Adds MSES z-scores from s5_se.
# ============================================================

message("\n--- s7.5: Multilevel model: Threshold ~ MAIA + MSES ---")

# Identify breath session per participant and filter
thresh_long_breath <- s5_thresh_long |>
  dplyr::left_join(
    s5_se |> dplyr::select(id, breath_ses, MSES_total_z, MSES_selfdoubt_z),
    by = "id"
  ) |>
  dplyr::filter(
    !is.na(Threshold),
    ses == breath_ses # keep breath session only
  )

cat(sprintf("Threshold long (breath sessions): %d rows, %d participants\n",
            nrow(thresh_long_breath),
            dplyr::n_distinct(thresh_long_breath$id)))

thresh_long_cc <- thresh_long_breath |>
  dplyr::filter(
    !is.na(Threshold), !is.na(MAIA_total_z),
    !is.na(MSES_total_z), !is.na(MSES_selfdoubt_z),
    !is.na(Salience), !is.na(Direction)
  )

cat(sprintf("Threshold complete cases: %d rows, %d participants\n",
            nrow(thresh_long_cc), dplyr::n_distinct(thresh_long_cc$id)))

# Base model
m_thresh_base <- lmerTest::lmer(
  Threshold ~ Salience + Direction + (1 | id),
  data = thresh_long_cc, REML = FALSE, control = ctrl
)

# MAIA main effect
m_thresh_maia <- lmerTest::lmer(
  Threshold ~ MAIA_total_z + Salience + Direction + (1 | id),
  data = thresh_long_cc, REML = FALSE, control = ctrl
)

# Add MSES controls
m_thresh_mses <- lmerTest::lmer(
  Threshold ~ MAIA_total_z + MSES_total_z + MSES_selfdoubt_z +
    Salience + Direction + (1 | id),
  data = thresh_long_cc, REML = FALSE, control = ctrl
)

cat("\nFixed effects (MAIA + MSES + Salience + Direction):\n")
print(round(summary(m_thresh_mses)$coefficients, 4))

lrt_thresh <- anova(m_thresh_base, m_thresh_maia, m_thresh_mses)
cat("\nLRT sequence (threshold models):\n")
print(lrt_thresh)

s7_thresh_mlm <- dplyr::bind_rows(
  .extract_mlm(m_thresh_maia, "MAIA + Salience + Direction", "breath"),
  .extract_mlm(m_thresh_mses, "MAIA + MSES + Salience + Direction", "breath")
)


# ============================================================
# Save outputs
# ============================================================

readr::write_csv(s7_conf_partials,
                 file.path(RESULTS_DIR, "s7_partial_correlations.csv"))
readr::write_csv(s7_thresh_partials,
                 file.path(RESULTS_DIR, "s7_threshold_partial_correlations.csv"))
readr::write_csv(s7_conf_mlm,
                 file.path(RESULTS_DIR, "s7_multilevel_confidence_results.csv"))
readr::write_csv(s7_thresh_mlm,
                 file.path(RESULTS_DIR, "s7_multilevel_threshold_results.csv"))

message("Saved: s7_partial_correlations.csv")
message("Saved: s7_threshold_partial_correlations.csv")
message("Saved: s7_multilevel_confidence_results.csv")
message("Saved: s7_multilevel_threshold_results.csv")
message("analysis_s7_maia_selfesteem.R complete.")
