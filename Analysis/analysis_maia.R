# ============================================================
# analysis_maia.R
# BCAT — MAIA: Sensibility vs. Sensitivity (H3A, H3B)
#
# Corresponds to: main text "MAIA Reflects Metacognitive Sensibility"
#
# NOTE: analysis_val_thresholds.R must be sourced before this file.
#       s{N}_maia, s{N}_thresh, s{N}_conf, s4_thresh_long, s5_thresh_long
#       are created there and used here.
#
# Assumes in environment (set by MainAnalysis.R + prior files):
#   s1s, s2s, s3s, s4s, s5s
#   s1l, s3l, s4l, s5_long_breath
#   s1_maia_A, s1_maia_B, s2_maia, s3_maia, s4_maia, s5_maia
#   s4_thresh, s5_thresh, s4_test, s5_test
#   s4_thresh_long, s5_thresh_long
#   s4_arousal, s5_arousal, meta_h4b
#   .extract_coef()    — defined in analysis_val_thresholds.R
#   .h4a_lrt()         — defined in analysis_arousal.R
#   .extract_h4c_coef() — defined in analysis_arousal.R
#   Helpers from utils.R: apply_bh_correction()
#   meta_analysis.R: run_h3_meta(), write_meta_results(), plot_meta_figures()
#   RESULTS_DIR, ANALYSIS_DIR
#
# Outputs:
#   table_maia.csv
#   table_maia_gating_moderation.csv
#   table_threshold_descriptives.csv
#   study4_bh_correction.csv
#   study5_bh_correction.csv
# ============================================================

source(file.path(ANALYSIS_DIR, "theme_bcat.R"))
message("\n========================================")
message("MAIA: SENSIBILITY VS. SENSITIVITY (H3A, H3B)")
message("========================================")


# ============================================================
# fit_maia_moderation()
#
# Exploratory: does MAIA moderate the Change -> Arousal relationship?
# Fits: Arousal ~ Change * MAIA_total_z + Change² + (1 | id)
#
# The interaction term (Change:MAIA_total_z) answers: do higher-MAIA
# individuals show stronger or weaker arousal responses to detected
# breathing changes?
#
# `data`: trial-level with MAIA_total_z merged in.
#         Required columns: id, Arousal, Change, MAIA_total_z.
# ============================================================
fit_maia_moderation <- function(data, study_label = "") {

  data <- data |>
    dplyr::filter(!is.na(Arousal), !is.na(Change), !is.na(MAIA_total_z)) |>
    dplyr::mutate(Change2 = Change^2)

  n_ppt <- dplyr::n_distinct(data$id)
  cat(sprintf("\n[%s] MAIA moderation: N=%d participants, %d trials\n",
              study_label, n_ppt, nrow(data)))

  if (n_ppt < 10) {
    message(sprintf("  [%s] Too few participants for MAIA moderation — skipping",
                    study_label))
    return(NULL)
  }

  ctrl <- lme4::lmerControl(optimizer = "bobyqa")

  m_base <- lmerTest::lmer(
    Arousal ~ Change + Change2 + (1 | id),
    data = data, REML = FALSE, control = ctrl)

  m_main <- lmerTest::lmer(
    Arousal ~ Change + Change2 + MAIA_total_z + (1 | id),
    data = data, REML = FALSE, control = ctrl)

  m_int <- lmerTest::lmer(
    Arousal ~ Change * MAIA_total_z + Change2 + (1 | id),
    data = data, REML = FALSE, control = ctrl)

  lrt <- anova(m_base, m_main, m_int)

  cf      <- summary(m_int)$coefficients
  int_row <- grep("Change:MAIA|MAIA.*Change", rownames(cf), value = TRUE)

  if (length(int_row) > 0) {
    b  <- cf[int_row[1], "Estimate"]
    se <- cf[int_row[1], "Std. Error"]
    t  <- cf[int_row[1], "t value"]
    df <- if ("df" %in% colnames(cf)) cf[int_row[1], "df"] else
            nobs(m_int) - length(lme4::fixef(m_int))
    p  <- cf[int_row[1], "Pr(>|t|)"]
    pr <- t / sqrt(t^2 + df)
    cat(sprintf(
      "  Change x MAIA_z: b=%.3f (SE=%.3f), partial r=%.3f, p=%.4f\n",
      b, se, pr, p))
  } else {
    message(sprintf("  [%s] Change:MAIA_total_z term not found", study_label))
  }

  list(base = m_base, main = m_main, int = m_int, lrt = lrt, study = study_label)
}


# ── H3A / H3B summary table ───────────────────────────────────

.bf01 <- function(maia_obj) {
  if (!is.null(maia_obj$H3B_BF) && !is.na(maia_obj$H3B_BF))
    round(1 / maia_obj$H3B_BF, 3)
  else NA_real_
}

.maia_val <- function(x) if (is.null(x)) NA_real_ else x

table_maia <- tibble::tibble(
  study        = c("Study1A", "Study2", "Study3", "Study4", "Study5"),
  n            = c(nrow(s1s), nrow(s2s), nrow(s3s), nrow(s4s), nrow(s5s)),
  H3A_r        = c(.maia_val(s1_maia_A$H3A_freq$estimate),
                   .maia_val(s2_maia$H3A_freq$estimate),
                   .maia_val(s3_maia$H3A_freq$estimate),
                   .maia_val(s4_maia$H3A_freq$estimate),
                   .maia_val(s5_maia$H3A_freq$estimate)),
  H3A_CI_lower = c(.maia_val(s1_maia_A$H3A_freq$conf.int[1]),
                   .maia_val(s2_maia$H3A_freq$conf.int[1]),
                   .maia_val(s3_maia$H3A_freq$conf.int[1]),
                   .maia_val(s4_maia$H3A_freq$conf.int[1]),
                   .maia_val(s5_maia$H3A_freq$conf.int[1])),
  H3A_CI_upper = c(.maia_val(s1_maia_A$H3A_freq$conf.int[2]),
                   .maia_val(s2_maia$H3A_freq$conf.int[2]),
                   .maia_val(s3_maia$H3A_freq$conf.int[2]),
                   .maia_val(s4_maia$H3A_freq$conf.int[2]),
                   .maia_val(s5_maia$H3A_freq$conf.int[2])),
  H3A_df       = c(.maia_val(s1_maia_A$H3A_freq$parameter),
                   .maia_val(s2_maia$H3A_freq$parameter),
                   .maia_val(s3_maia$H3A_freq$parameter),
                   .maia_val(s4_maia$H3A_freq$parameter),
                   .maia_val(s5_maia$H3A_freq$parameter)),
  H3A_p        = c(.maia_val(s1_maia_A$H3A_freq$p.value),
                   .maia_val(s2_maia$H3A_freq$p.value),
                   .maia_val(s3_maia$H3A_freq$p.value),
                   .maia_val(s4_maia$H3A_freq$p.value),
                   .maia_val(s5_maia$H3A_freq$p.value)),
  H3B_r        = c(.maia_val(s1_maia_A$H3B_freq$estimate),
                   .maia_val(s2_maia$H3B_freq$estimate),
                   .maia_val(s3_maia$H3B_freq$estimate),
                   .maia_val(s4_maia$H3B_freq$estimate),
                   .maia_val(s5_maia$H3B_freq$estimate)),
  H3B_CI_lower = c(.maia_val(s1_maia_A$H3B_freq$conf.int[1]),
                   .maia_val(s2_maia$H3B_freq$conf.int[1]),
                   .maia_val(s3_maia$H3B_freq$conf.int[1]),
                   .maia_val(s4_maia$H3B_freq$conf.int[1]),
                   .maia_val(s5_maia$H3B_freq$conf.int[1])),
  H3B_CI_upper = c(.maia_val(s1_maia_A$H3B_freq$conf.int[2]),
                   .maia_val(s2_maia$H3B_freq$conf.int[2]),
                   .maia_val(s3_maia$H3B_freq$conf.int[2]),
                   .maia_val(s4_maia$H3B_freq$conf.int[2]),
                   .maia_val(s5_maia$H3B_freq$conf.int[2])),
  H3B_df       = c(.maia_val(s1_maia_A$H3B_freq$parameter),
                   .maia_val(s2_maia$H3B_freq$parameter),
                   .maia_val(s3_maia$H3B_freq$parameter),
                   .maia_val(s4_maia$H3B_freq$parameter),
                   .maia_val(s5_maia$H3B_freq$parameter)),
  H3B_p        = c(.maia_val(s1_maia_A$H3B_freq$p.value),
                   .maia_val(s2_maia$H3B_freq$p.value),
                   .maia_val(s3_maia$H3B_freq$p.value),
                   .maia_val(s4_maia$H3B_freq$p.value),
                   .maia_val(s5_maia$H3B_freq$p.value)),
  H3B_BF01     = c(.bf01(s1_maia_A), .bf01(s2_maia), .bf01(s3_maia),
                   .bf01(s4_maia),   .bf01(s5_maia))
)
readr::write_csv(table_maia, file.path(RESULTS_DIR, "table_maia.csv"))
message("Saved: table_maia.csv")


# ── MAIA meta-analysis (H3A, H3B) ─────────────────────────────

maia_table_meta <- table_maia |>
  dplyr::mutate(
    n_Confidence      = n,
    # Study 3 excluded from pooled H3B: uses mean accuracy as proxy,
    # not a staircase-derived threshold — not comparable across studies.
    n_Threshold       = dplyr::if_else(study == "Study3", NA_integer_, as.integer(n)),
    r_Confidence_MAIA = H3A_r,
    r_Threshold_MAIA  = dplyr::if_else(study == "Study3", NA_real_, H3B_r)) |>
  dplyr::select(-n)
meta_h3 <- run_h3_meta(maia_table_meta)


# ── MAIA moderation (exploratory; Studies 3, 4, 5) ────────────

s3_maia_mod <- fit_maia_moderation(
  s3l |> dplyr::left_join(
    s3s |> dplyr::select(id, MAIA_total_z), by = "id"),
  study_label = "Study3")
s4_maia_mod <- fit_maia_moderation(
  s4l |> dplyr::left_join(
    s4s |> dplyr::select(id, MAIA_total_z), by = "id"),
  study_label = "Study4")
s5_maia_mod <- fit_maia_moderation(
  s5_long_breath |> dplyr::left_join(
    s5s |> dplyr::select(id, MAIA_total_z), by = "id"),
  study_label = "Study5")


# ============================================================
# fit_maia_gating_moderation()
#
# Exploratory: does MAIA moderate the awareness-gating effect?
# Tests the three-way Change × Accuracy × MAIA_z interaction.
#
# Two theoretically distinct questions:
#   (a) Accuracy × MAIA_z: do high-MAIA individuals show larger
#       hit/miss differences in baseline arousal (intercept gating)?
#   (b) Change × Accuracy × MAIA_z: do high-MAIA individuals show
#       stronger magnitude-scaling on detected vs missed trials?
#
# LRT chain:
#   m_gate    : Arousal ~ Change * Accuracy + Change2 (standard gating)
#   m_gate_m  : + MAIA_total_z                       (MAIA as covariate)
#   m_gate_ax : + Accuracy * MAIA_total_z            (gating × MAIA)
#   m_gate_3w : + Change * Accuracy * MAIA_total_z   (full three-way)
#
# `data`: trial-level; columns id, Arousal, Change, Accuracy, MAIA_total_z.
# ============================================================
fit_maia_gating_moderation <- function(data, study_label = "") {

  data <- data |>
    dplyr::filter(
      !is.na(Arousal), !is.na(Change),
      !is.na(Accuracy), !is.na(MAIA_total_z)
    ) |>
    dplyr::mutate(
      Change2  = Change^2,
      Accuracy = as.numeric(Accuracy)   # 0/1 numeric for interaction
    )

  n_ppt <- dplyr::n_distinct(data$id)
  cat(sprintf(
    "\n[%s] MAIA gating moderation: N=%d participants, %d trials\n",
    study_label, n_ppt, nrow(data)
  ))

  if (n_ppt < 10) {
    message(sprintf("  [%s] Too few participants — skipping", study_label))
    return(NULL)
  }

  ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

  m_gate <- lmerTest::lmer(
    Arousal ~ Change * Accuracy + Change2 + (1 | id),
    data = data, REML = FALSE, control = ctrl)

  m_gate_m <- lmerTest::lmer(
    Arousal ~ Change * Accuracy + MAIA_total_z + Change2 + (1 | id),
    data = data, REML = FALSE, control = ctrl)

  m_gate_ax <- lmerTest::lmer(
    Arousal ~ Change * Accuracy + Accuracy * MAIA_total_z + Change2 + (1 | id),
    data = data, REML = FALSE, control = ctrl)

  m_gate_3w <- lmerTest::lmer(
    Arousal ~ Change * Accuracy * MAIA_total_z + Change2 + (1 | id),
    data = data, REML = FALSE, control = ctrl)

  lrt <- anova(m_gate, m_gate_m, m_gate_ax, m_gate_3w)

  cf <- summary(m_gate_3w)$coefficients

  # Helper: extract named row safely using fixed string matching
  # lmerTest names interaction terms with colons: "A:B", "A:B:C"
  .get <- function(term_name) {
    row <- rownames(cf)[rownames(cf) == term_name]
    if (length(row) == 0) {
      # fallback: partial match without regex
      row <- rownames(cf)[startsWith(rownames(cf), term_name)]
    }
    if (length(row) == 0) return(c(b = NA_real_, se = NA_real_,
                                   t = NA_real_,  p = NA_real_))
    x <- cf[row[1], c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
    c(b = x[1], se = x[2], t = x[3], p = x[4])
  }

  acc_maia        <- .get("Accuracy:MAIA_total_z")
  change_acc_maia <- .get("Change:Accuracy:MAIA_total_z")

  cat(sprintf(
    "  Accuracy × MAIA_z: b=%.3f (SE=%.3f), t=%.2f, p=%.4f\n",
    acc_maia["b"], acc_maia["se"], acc_maia["t"], acc_maia["p"]
  ))
  cat(sprintf(
    "  Change × Accuracy × MAIA_z: b=%.3f (SE=%.3f), t=%.2f, p=%.4f\n",
    change_acc_maia["b"], change_acc_maia["se"],
    change_acc_maia["t"], change_acc_maia["p"]
  ))
  cat("  LRT (gating → +MAIA → +Acc×MAIA → +3-way):\n")
  print(lrt)

  list(
    m_gate    = m_gate,
    m_gate_m  = m_gate_m,
    m_gate_ax = m_gate_ax,
    m_gate_3w = m_gate_3w,
    lrt       = lrt,
    study     = study_label
  )
}

# ── Run for Studies 4 and 5 (have trial-level Accuracy + MAIA) ─
# Study 3 excluded: uses mean accuracy as proxy, not trial-level.

s4_maia_gating <- fit_maia_gating_moderation(
  s4l |> dplyr::left_join(
    s4s |> dplyr::select(id, MAIA_total_z), by = "id"),
  study_label = "Study4")

s5_maia_gating <- fit_maia_gating_moderation(
  s5_long_breath |> dplyr::left_join(
    s5s |> dplyr::select(id, MAIA_total_z), by = "id"),
  study_label = "Study5")

# ── Save key terms to CSV ─────────────────────────────────────
.extract_gating_terms <- function(obj) {
  if (is.null(obj)) return(NULL)
  cf  <- summary(obj$m_gate_3w)$coefficients
  lrt <- obj$lrt

  tibble::tibble(
    study = obj$study,
    term  = rownames(cf),
    b     = cf[, "Estimate"],
    se    = cf[, "Std. Error"],
    t     = cf[, "t value"],
    df    = if ("df" %in% colnames(cf)) cf[, "df"] else NA_real_,
    p     = cf[, "Pr(>|t|)"],
    partial_r = t / sqrt(t^2 + df)
  ) |>
    dplyr::mutate(
      lrt_acc_maia_p  = lrt[["Pr(>Chisq)"]][3],   # +Acc×MAIA step
      lrt_threeway_p  = lrt[["Pr(>Chisq)"]][4],   # +Change×Acc×MAIA step
      dplyr::across(where(is.numeric), \(x) round(x, 4))
    )
}

table_maia_gating <- dplyr::bind_rows(
  .extract_gating_terms(s4_maia_gating),
  .extract_gating_terms(s5_maia_gating)
)

readr::write_csv(
  table_maia_gating,
  file.path(RESULTS_DIR, "table_maia_gating_moderation.csv")
)
message("Saved: table_maia_gating_moderation.csv")

table_threshold_descriptives <- dplyr::bind_rows(
  s4_thresh_long |>
    dplyr::filter(!is.na(Threshold)) |>
    dplyr::group_by(Salience, Direction) |>
    dplyr::summarise(
      M_threshold  = mean(Threshold, na.rm = TRUE),
      SD_threshold = sd(Threshold,   na.rm = TRUE),
      n_ppts       = dplyr::n_distinct(id), .groups = "drop") |>
    dplyr::mutate(study = "Study4"),
  s5_thresh_long |>
    dplyr::filter(!is.na(Threshold)) |>
    dplyr::group_by(ses, Salience, Direction) |>
    dplyr::summarise(
      M_threshold  = mean(Threshold, na.rm = TRUE),
      SD_threshold = sd(Threshold,   na.rm = TRUE),
      n_ppts       = dplyr::n_distinct(id), .groups = "drop") |>
    dplyr::mutate(study = "Study5")
) |> dplyr::select(study, dplyr::everything())
readr::write_csv(table_threshold_descriptives,
                 file.path(RESULTS_DIR, "table_threshold_descriptives.csv"))
message("Saved: table_threshold_descriptives.csv")


# ── BH correction (Studies 4 and 5; pre-registered only) ──────

.h4b_p <- function(arousal_obj) {
  cf       <- summary(arousal_obj$H4B_final)$coefficients
  int_term <- grep("Change.*Accuracy|Accuracy.*Change", rownames(cf), value = TRUE)
  cf[int_term[1], "Pr(>|t|)"]
}

s4_pvals <- c(
  H1  = .extract_coef(s4_thresh$H1,   "DirectionSlower")[["p"]],
  H2  = .extract_coef(s4_thresh$H2,   "SalienceHigh")[["p"]],
  H3A = .maia_val(s4_maia$H3A_freq$p.value),
  H3B = .maia_val(s4_maia$H3B_freq$p.value),
  H4A = .h4a_lrt(s4_arousal$H4A_null, s4_arousal$H4A_quad)[["p"]],
  H5  = .extract_coef(s4_test$H5_sal, "SalienceHigh", p_col = "Pr(>|z|)")[["p"]],
  H4B = .h4b_p(s4_arousal),
  H4C = .extract_h4c_coef(s4_arousal$H4C_final, "Study4")[["p"]]
)
s4_bh <- apply_bh_correction(s4_pvals, study_label = "Study4")
readr::write_csv(s4_bh, file.path(RESULTS_DIR, "study4_bh_correction.csv"))
message("Saved: study4_bh_correction.csv")

s5_pvals <- c(
  H1  = .extract_coef(s5_thresh$H1,   "DirectionSlower")[["p"]],
  H2  = .extract_coef(s5_thresh$H2,   "SalienceHigh")[["p"]],
  H3A = .maia_val(s5_maia$H3A_freq$p.value),
  H3B = .maia_val(s5_maia$H3B_freq$p.value),
  H4A = .h4a_lrt(s5_arousal$H4A_null, s5_arousal$H4A_quad)[["p"]],
  H5  = .extract_coef(s5_test$H5_sal, "SalienceHigh", p_col = "Pr(>|z|)")[["p"]],
  H4B = .h4b_p(s5_arousal),
  H4C = .extract_h4c_coef(s5_arousal$H4C_final, "Study5")[["p"]]
)
s5_bh <- apply_bh_correction(s5_pvals, study_label = "Study5")
readr::write_csv(s5_bh, file.path(RESULTS_DIR, "study5_bh_correction.csv"))
message("Saved: study5_bh_correction.csv")


# ── Write meta results ────────────────────────────────────────
write_meta_results(meta_h4b, meta_h3)
plot_meta_figures(meta_h4b, meta_h3)
