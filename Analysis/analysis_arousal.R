# ============================================================
# analysis_arousal.R
# BCAT — Awareness Gates Arousal Transfer (H4A, H4B, H4C)
#
# Corresponds to: main text "Awareness Gates Arousal Transfer" section
#
# Assumes in environment (set by MainAnalysis.R):
#   s1l, s1s, s2l, s2s, s3l, s3s, s4l, s4s, s5l, s5s
#   s5_long_breath
#   RESULTS_DIR, ANALYSIS_DIR
#   Helpers from utils.R: extract_h4b(), partial_r_from_t(), fmt_lmer()
#   meta_analysis.R: run_h4b_meta(), write_meta_results(), plot_meta_figures()
#
# Creates in environment (used downstream by analysis_maia.R):
#   s1_arousal_A, s2_arousal, s3_arousal, s4_arousal, s5_arousal
#   .h4a_lrt(), .extract_h4c_coef() — local helpers also used by analysis_maia.R
#   meta_h4b
#
# Outputs:
#   table_arousal.csv
#   meta_h4b_sensitivity_noS3.csv
#   table_h4c_exploratory.csv
# ============================================================

message("\n========================================")
message("AWARENESS GATES AROUSAL TRANSFER (H4A, H4B, H4C)")
message("========================================")


# ============================================================
# fit_arousal_models()
#
# Core arousal LMM family (H4A, H4B, H4C).
# Follows the unified pre-registered structure:
#   H4A: Change main effect (linear + quadratic)
#   H4B: Change x Accuracy moderation (awareness gating)
#   H4C: Change x Condition/Group moderation (instruction effect)
#
# Model selection via LRT (ML); final models refitted with REML.
# R² via MuMIn::r.squaredGLMM().
#
# `data`:        trial-level; canonical columns id, Change, Arousal, Accuracy.
# `group_var`:   column name for H4C grouping (NULL to skip).
# `quadratic`:   include Change² (default TRUE — consistent with preregistration).
# `study_label`: string for console output.
# ============================================================
fit_arousal_models <- function(data,
                               group_var   = NULL,
                               quadratic   = TRUE,
                               study_label = "") {

  data <- data |>
    dplyr::filter(!is.na(Arousal), !is.na(Change)) |>
    dplyr::mutate(
      Change2  = Change^2,
      Accuracy = as.numeric(Accuracy)
    )

  cat(sprintf("\n[%s] Arousal models: N=%d trials, %d participants\n",
              study_label, nrow(data), dplyr::n_distinct(data$id)))

  # ── H4A: main effect of Change ───────────────────────────────
  m0  <- lmerTest::lmer(Arousal ~ 1      + (1 | id), data = data, REML = FALSE)
  mA  <- lmerTest::lmer(Arousal ~ Change + (1 | id), data = data, REML = FALSE)

  if (quadratic) {
    mAq <- lmerTest::lmer(Arousal ~ Change + Change2 + (1 | id),
                       data = data, REML = FALSE)
  } else {
    mAq <- mA
  }

  lrt_H4A <- anova(m0, mA, mAq)

  # ── H4B: Change × Accuracy ───────────────────────────────────
  data_b <- data |> dplyr::filter(!is.na(Accuracy))

  if (quadratic) {
    mB <- lmerTest::lmer(Arousal ~ Change * Accuracy + Change2 + (1 | id),
                      data = data_b, REML = FALSE)
  } else {
    mB <- lmerTest::lmer(Arousal ~ Change * Accuracy + (1 | id),
                      data = data_b, REML = FALSE)
  }
  lrt_H4B <- anova(update(mAq, data = data_b, REML = FALSE), mB)

  results <- list(
    H4A_null   = m0,
    H4A_linear = mA,
    H4A_quad   = mAq,
    H4B        = mB,
    lrt_H4A    = lrt_H4A,
    lrt_H4B    = lrt_H4B,
    study      = study_label
  )

  cat(sprintf("  H4A Change: %s\n", fmt_lmer(mAq, "Change")))
  cat(sprintf("  H4B Change x Accuracy: %s\n",
              fmt_lmer(mB, grep("Change.*Accuracy|Accuracy.*Change",
                                  rownames(summary(mB)$coefficients),
                                  value = TRUE)[1])))

  # ── H4C: Change × Group/Condition ────────────────────────────
  if (!is.null(group_var) && group_var %in% names(data)) {
    data[[".group"]] <- factor(data[[group_var]])

    if (quadratic) {
      mC <- lmerTest::lmer(Arousal ~ Change * .group + Change2 + (1 | id),
                        data = data, REML = FALSE)
    } else {
      mC <- lmerTest::lmer(Arousal ~ Change * .group + (1 | id),
                        data = data, REML = FALSE)
    }

    results$H4C     <- mC
    results$lrt_H4C <- anova(mAq, mC)

    int_term <- grep("\\.group", rownames(summary(mC)$coefficients), value = TRUE)
    int_term <- int_term[grep(":", int_term)]
    if (length(int_term) > 0)
      cat(sprintf("  H4C Change x %s: %s\n", group_var, fmt_lmer(mC, int_term[1])))
  }

  # ── REML refits for final reporting ──────────────────────────
  results$H4A_final <- update(mAq, REML = TRUE)
  results$H4B_final <- update(mB,  REML = TRUE)
  if (!is.null(results$H4C))
    results$H4C_final <- update(results$H4C, REML = TRUE)

  results$R2_H4A <- tryCatch(
    MuMIn::r.squaredGLMM(mAq), error = function(e) c(R2m = NA, R2c = NA))
  results$R2_H4B <- tryCatch(
    MuMIn::r.squaredGLMM(mB),  error = function(e) c(R2m = NA, R2c = NA))

  cat(sprintf("  R2m H4A=%.3f  H4B=%.3f\n",
              results$R2_H4A[1], results$R2_H4B[1]))

  results
}


# ── Run arousal models across all studies ─────────────────────

s1_arousal_A <- fit_arousal_models(
  dplyr::filter(s1l, Group == "TaskA"), study_label = "Study1_TaskA")
s2_arousal   <- fit_arousal_models(s2l, study_label = "Study2")
s3_arousal   <- fit_arousal_models(
  dplyr::filter(s3l, !is.na(Accuracy)),
  group_var = "Run", study_label = "Study3")
s4_arousal   <- fit_arousal_models(
  s4l, group_var = "Group", study_label = "Study4")
s5_arousal   <- fit_arousal_models(s5_long_breath, study_label = "Study5")


# ── H4C: interoceptive specificity (Studies 4 and 5) ─────────
#
# Primary: between-groups ses1 only.
# Exploratory: within-person visual-first participants
# (note: Condition and ses are collinear in that subset).

paced_conds <- c("highSalienceAcc", "highSalienceDec",
                 "lowSalienceAcc",  "lowSalienceDec")

visual_first_ids <- s5l |>
  dplyr::filter(ses == "ses1", Condition == "visual") |>
  dplyr::distinct(id) |>
  dplyr::pull(id)

s5_long_h4c_s1 <- s5l |>
  dplyr::filter(ses == "ses1",
                taskCondition %in% paced_conds,
                !is.na(Arousal), !is.na(Change), !is.na(Condition)) |>
  dplyr::mutate(Change2   = Change^2,
                Condition = factor(Condition, levels = c("breath", "visual")))

m_H4C <- lmerTest::lmer(
  Arousal ~ Change * Condition + Change2 + (1 | id),
  data = s5_long_h4c_s1, REML = FALSE)

s5_long_h4c_vf <- s5l |>
  dplyr::filter(id %in% visual_first_ids,
                taskCondition %in% paced_conds,
                !is.na(Arousal), !is.na(Change), !is.na(Condition)) |>
  dplyr::mutate(Change2   = Change^2,
                Condition = factor(Condition, levels = c("breath", "visual")))

m_H4C_exp <- lmerTest::lmer(
  Arousal ~ Change * Condition + Change2 + (1 | id),
  data = s5_long_h4c_vf, REML = FALSE)

s5_arousal$H4C           <- m_H4C
s5_arousal$H4C_final     <- update(m_H4C,     REML = TRUE)
s5_arousal$H4C_exp       <- m_H4C_exp
s5_arousal$H4C_exp_final <- update(m_H4C_exp, REML = TRUE)


# ── Cross-study arousal summary table ─────────────────────────

# Returns key stats from the H4C interaction term
.extract_h4c_coef <- function(model, study_label = "") {
  na_out <- c(b = NA_real_, se = NA_real_, t = NA_real_,
              df = NA_real_, p = NA_real_, partial_r = NA_real_)
  if (is.null(model)) return(na_out)
  cf <- tryCatch(summary(model)$coefficients, error = function(e) NULL)
  if (is.null(cf)) return(na_out)
  row_nm <- grep(
    "Change.*[Cc]ondition|Change.*[Gg]roup|[Cc]ondition.*Change|[Gg]roup.*Change",
    rownames(cf), value = TRUE)[1]
  if (is.na(row_nm)) return(na_out)
  b   <- cf[row_nm, "Estimate"]
  se  <- cf[row_nm, "Std. Error"]
  t   <- cf[row_nm, "t value"]
  df  <- if ("df" %in% colnames(cf)) cf[row_nm, "df"] else
           nobs(model) - length(lme4::fixef(model))
  p   <- cf[row_nm, if ("Pr(>|t|)" %in% colnames(cf)) "Pr(>|t|)" else "Pr(>|z|)"]
  c(b = b, se = se, t = t, df = df, p = p, partial_r = t / sqrt(t^2 + df))
}

# Returns c(chi2, df, p) from a null vs. quadratic LRT
.h4a_lrt <- function(null_mod, quad_mod) {
  lrt <- tryCatch(anova(null_mod, quad_mod), error = function(e) NULL)
  if (is.null(lrt))
    return(c(chi2 = NA_real_, df = NA_real_, p = NA_real_))
  c(chi2 = lrt[["Chisq"]][2], df = lrt[["Df"]][2], p = lrt[["Pr(>Chisq)"]][2])
}

h4b_table <- dplyr::bind_rows(
  as.data.frame(t(extract_h4b(s1_arousal_A$H4B_final, "Study1A"))) |>
    dplyr::mutate(study = "Study1A", n = nrow(s1s)),
  as.data.frame(t(extract_h4b(s2_arousal$H4B_final,   "Study2"))) |>
    dplyr::mutate(study = "Study2",  n = nrow(s2s)),
  as.data.frame(t(extract_h4b(s3_arousal$H4B_final,   "Study3"))) |>
    dplyr::mutate(study = "Study3",  n = nrow(s3s)),
  as.data.frame(t(extract_h4b(s4_arousal$H4B_final,   "Study4"))) |>
    dplyr::mutate(study = "Study4",  n = nrow(s4s)),
  as.data.frame(t(extract_h4b(s5_arousal$H4B_final,   "Study5"))) |>
    dplyr::mutate(study = "Study5",  n = nrow(s5s))
)

.h4c_s4 <- .extract_h4c_coef(s4_arousal$H4C_final, "Study4")
.h4c_s5 <- .extract_h4c_coef(s5_arousal$H4C_final, "Study5")

.h4a <- list(
  s1 = .h4a_lrt(s1_arousal_A$H4A_null, s1_arousal_A$H4A_quad),
  s2 = .h4a_lrt(s2_arousal$H4A_null,   s2_arousal$H4A_quad),
  s3 = .h4a_lrt(s3_arousal$H4A_null,   s3_arousal$H4A_quad),
  s4 = .h4a_lrt(s4_arousal$H4A_null,   s4_arousal$H4A_quad),
  s5 = .h4a_lrt(s5_arousal$H4A_null,   s5_arousal$H4A_quad)
)

table_arousal <- tibble::tibble(
  study         = c("Study1A", "Study2", "Study3", "Study4", "Study5"),
  n             = h4b_table$n,
  H4A_chi2      = c(.h4a$s1["chi2"], .h4a$s2["chi2"], .h4a$s3["chi2"],
                    .h4a$s4["chi2"], .h4a$s5["chi2"]),
  H4A_df        = c(.h4a$s1["df"],   .h4a$s2["df"],   .h4a$s3["df"],
                    .h4a$s4["df"],   .h4a$s5["df"]),
  H4A_p_LRT     = c(.h4a$s1["p"],    .h4a$s2["p"],    .h4a$s3["p"],
                    .h4a$s4["p"],    .h4a$s5["p"]),
  H4B_b         = h4b_table$b,
  H4B_SE        = h4b_table$se,
  H4B_t         = h4b_table$t,
  H4B_df        = h4b_table$df,
  H4B_p         = h4b_table$p,
  H4B_partial_r = h4b_table$partial_r,
  H4C_b         = c(NA_real_, NA_real_, NA_real_, .h4c_s4["b"],        .h4c_s5["b"]),
  H4C_SE        = c(NA_real_, NA_real_, NA_real_, .h4c_s4["se"],       .h4c_s5["se"]),
  H4C_t         = c(NA_real_, NA_real_, NA_real_, .h4c_s4["t"],        .h4c_s5["t"]),
  H4C_df        = c(NA_real_, NA_real_, NA_real_, .h4c_s4["df"],       .h4c_s5["df"]),
  H4C_p         = c(NA_real_, NA_real_, NA_real_, .h4c_s4["p"],        .h4c_s5["p"]),
  H4C_partial_r = c(NA_real_, NA_real_, NA_real_, .h4c_s4["partial_r"], .h4c_s5["partial_r"])
)

readr::write_csv(table_arousal, file.path(RESULTS_DIR, "table_arousal.csv"))
message("Saved: table_arousal.csv")

readr::write_csv(
  tibble::tibble(
    study  = "Study5_VisualFirst",
    H4C_b  = .extract_h4c_coef(s5_arousal$H4C_exp_final, "Study5_exp")["b"],
    H4C_SE = .extract_h4c_coef(s5_arousal$H4C_exp_final, "Study5_exp")["se"],
    H4C_p  = .extract_h4c_coef(s5_arousal$H4C_exp_final, "Study5_exp")["p"],
    note   = "Visual-first participants only; within-person; exploratory"),
  file.path(RESULTS_DIR, "table_h4c_exploratory.csv"))
message("Saved: table_h4c_exploratory.csv")


# ── Meta-analysis (H4B, k=5) ──────────────────────────────────

study_results <- list(
  Study1_TaskA = s1_arousal_A,
  Study2       = s2_arousal,
  Study3       = s3_arousal,
  Study4       = s4_arousal,
  Study5       = s5_arousal
)

study_ns <- list(
  Study1_TaskA = 181L,
  Study2 = 166L, Study3 = 103L, Study4 = 126L, Study5 = 206L
)

meta_h4b <- run_h4b_meta(study_results, study_ns)

# Sensitivity: exclude Study 3 (fixed large magnitudes, SAM scale)
meta_h4b_noS3 <- run_h4b_meta(
  study_results[names(study_results) != "Study3"],
  study_ns[names(study_ns) != "Study3"])
if (!is.null(meta_h4b_noS3)) {
  rma_s <- meta_h4b_noS3$rma
  readr::write_csv(
    tibble::tibble(
      analysis  = "H4B sensitivity: Study3 excluded",
      k_studies = rma_s$k,
      r_pooled  = tanh(rma_s$beta[[1]]),
      r_lower   = tanh(rma_s$ci.lb),
      r_upper   = tanh(rma_s$ci.ub),
      p_value   = rma_s$pval,
      I2_pct    = rma_s$I2,
      tau2      = rma_s$tau2,
      Q_stat    = rma_s$QE,
      Q_p       = rma_s$QEp),
    file.path(RESULTS_DIR, "meta_h4b_sensitivity_noS3.csv"))
  message("Saved: meta_h4b_sensitivity_noS3.csv")
}
