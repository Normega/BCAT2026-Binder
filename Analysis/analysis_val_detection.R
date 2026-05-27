# ============================================================
# analysis_val_detection.R
# BCAT — Task Validation: Change² Dose-Response (all studies)
#
# Corresponds to: main text "Task Validation" paragraph (Table 2)
#
# Assumes in environment (set by MainAnalysis.R):
#   s1l, s1s, s2l, s2s, s3l, s3s, s4l, s4s, s5l, s5s
#   RESULTS_DIR, ANALYSIS_DIR
#   Helpers from utils.R: extract_change2_results(), add_odds_ratios()
#
# Creates in environment (used downstream):
#   s1_det_A, s1_det_B — detection model objects for Study 1 Tasks A/B
#   s2_det, s3_det, s4_det, s5_det — detection model objects per study
#   table_detection_change2 — cross-study Change² summary tibble
#
# NOTE: table_validation.csv (assembled in analysis_val_thresholds.R)
#       depends on the s{N}_det objects created here. Run this file first.
#
# Output: table_detection_change2.csv
# ============================================================

message("\n========================================")
message("TASK VALIDATION: Change² detection models")
message("========================================")


# ============================================================
# fit_detection_models()
#
# Hierarchical detection GLMM matching Kyle's dissertation Table 2.1.
# Tests whether Change (and Change²) predicts trial-level detection
# accuracy. Salience models added when Salience is available.
#
# `data`: trial-level; canonical columns id, Accuracy, Change.
#   Salience column optional — models added if present.
# `study_label`: string for console output.
#
# Returns named list of glmer objects and LRT table.
# ============================================================
fit_detection_models <- function(data, study_label = "") {

  data <- data |>
    dplyr::filter(!is.na(Accuracy), !is.na(Change)) |>
    dplyr::mutate(
      Accuracy = as.integer(Accuracy),
      Change2  = Change^2
    )

  cat(sprintf("\n[%s] Detection models: N=%d trials, %d participants\n",
              study_label, nrow(data), dplyr::n_distinct(data$id)))

  # bobyqa is more robust than the default optimizer for binary glmer
  # with small per-person trial counts (e.g. Study 1: 15 trials/participant).
  # Singular fits (random intercept variance -> 0) are expected in this
  # situation — fixed effects are still valid; flag but do not abort.
  glmer_ctrl <- lme4::glmerControl(optimizer = "bobyqa",
                                    optCtrl   = list(maxfun = 2e5))

  .fit_glmer <- function(formula, data, ctrl) {
    m <- lme4::glmer(formula, data = data,
                     family  = binomial,
                     control = ctrl)
    if (lme4::isSingular(m)) {
      message(sprintf(
        "  [%s] Singular fit — random intercept variance near zero. ",
        study_label),
        "Fixed effects are unaffected; interpret random effects cautiously.")
    }
    m
  }

  m0  <- .fit_glmer(Accuracy ~ 1               + (1 | id), data, glmer_ctrl)
  m1  <- .fit_glmer(Accuracy ~ Change           + (1 | id), data, glmer_ctrl)
  m2  <- .fit_glmer(Accuracy ~ Change + Change2 + (1 | id), data, glmer_ctrl)

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

    m3 <- .fit_glmer(
      Accuracy ~ Change + Change2 + Salience + (1 | id), data, glmer_ctrl)
    m4 <- .fit_glmer(
      Accuracy ~ (Change + Change2) * Salience + (1 | id), data, glmer_ctrl)
    m5 <- tryCatch(
      .fit_glmer(
        Accuracy ~ (Change + Change2) * Salience + (1 + Salience | id),
        data, glmer_ctrl),
      error = function(e) {
        message("  Random-slope model failed; using fixed-effect interaction")
        m4
      }
    )
    results$salience_add <- m3
    results$salience_int <- m4
    results$salience_rs  <- m5
    results$lrt_salience <- anova(m2, m3, m4)
  }

  results
}


# ── Run detection models across all studies ───────────────────

s1_det_A <- fit_detection_models(
  dplyr::filter(s1l, Group == "TaskA"), study_label = "Study1_TaskA")
s1_det_B <- fit_detection_models(
  dplyr::filter(s1l, Group == "TaskB"), study_label = "Study1_TaskB")
s2_det   <- fit_detection_models(s2l, study_label = "Study2")
s3_det   <- fit_detection_models(s3l, study_label = "Study3")
s4_det   <- fit_detection_models(s4l, study_label = "Study4")
s5_det   <- fit_detection_models(s5l, study_label = "Study5")


# ── Extract Change² betas and LRT statistics ──────────────────

table_detection_change2 <- dplyr::bind_rows(
  extract_change2_results(s1_det_A, "Study1_TaskA", n_participants = nrow(s1s)),
  extract_change2_results(s1_det_B, "Study1_TaskB", n_participants = nrow(s1s)),
  extract_change2_results(s2_det,   "Study2",       n_participants = nrow(s2s)),
  extract_change2_results(s3_det,   "Study3",       n_participants = nrow(s3s)),
  extract_change2_results(s4_det,   "Study4",       n_participants = nrow(s4s)),
  extract_change2_results(s5_det,   "Study5",       n_participants = nrow(s5s))
)

# Add odds ratios (OR = exp(b); 95% CI = exp(b +/- 1.96*SE))
# NOTE: Change is signed in Studies 4/5 (negative = faster breathing);
# OR_Change < 1 there reflects directional coding, not a reversal.
# OR_Change2 is the primary psychophysical validity indicator.
table_detection_change2 <- add_odds_ratios(
  table_detection_change2, terms = c("Change", "Change2"))

readr::write_csv(table_detection_change2,
                 file.path(RESULTS_DIR, "table_detection_change2.csv"))
message("Saved: table_detection_change2.csv")
