# ============================================================
# belt_salience_followup.R
# Belt salience follow-up analyses (exploratory).
# Sources after analysis_belt.R; uses qcFull and s5_long_breath.
# Outputs: table_belt_salience_followup.csv
# ============================================================



# Paths — adjust to your environment


# ── 1. Filter qcFile to analysis sample ──────────────────────
# qcFull loaded by analysis_belt.R (FullResults sheet).
# Columns used: condition, correct_breaths, lag_flag, dur_b1-b4,
#   dur_b1_vs_b4, belt_quality, delta, salience, id, run.

qc <- qcFull |>
  dplyr::filter(
    condition       == "breath",
    correct_breaths == 1,
    lag_flag        == 0,
    !is.na(dur_b1),
    dplyr::between(dur_b1, 2, 7),
    belt_quality    != "unusable"
  ) |>
  dplyr::mutate(
    abs_delta        = abs(delta),
    abs_dur_b1_vs_b4 = abs(dur_b1_vs_b4),
    early_change     = dur_b2 - dur_b1,
    late_change      = dur_b4 - dur_b3,
    mid_change       = dur_b3 - dur_b2,
    # OLS slope of dur across 4 breaths (x = 1:4, mean_x = 2.5)
    linear_slope     = (
      (1 - 2.5) * dur_b1 + (2 - 2.5) * dur_b2 +
      (3 - 2.5) * dur_b3 + (4 - 2.5) * dur_b4
    ) / 5,
    onset_ratio      = dplyr::if_else(
      abs(dur_b1_vs_b4) > 0.01,
      early_change / dur_b1_vs_b4,
      NA_real_),
    salience_num     = dplyr::if_else(salience == "High", 1L, 0L)
  )

cat("R2 breath trials:", nrow(qc), "\n")
cat("Participants:",     dplyr::n_distinct(qc$id), "\n\n")
# ── 2. Effect size for the original Salience model ───────────────────────────
# Replicate original model, then compute partial r and partial eta^2

m_base <- lmerTest::lmer(
  abs_dur_b1_vs_b4 ~ abs_delta + (1 | id),
  data = qc, REML = FALSE
)
m_sal <- lmerTest::lmer(
  abs_dur_b1_vs_b4 ~ salience + abs_delta + (1 | id),
  data = qc, REML = FALSE
)

sal_coef  <- summary(m_sal)$coefficients["salienceLow", ]
t_sal     <- sal_coef["t value"]
df_resid  <- sal_coef["df"]
r_partial <- t_sal / sqrt(t_sal^2 + df_resid)
eta2_part <- t_sal^2 / (t_sal^2 + df_resid)   # partial eta^2

r2_base   <- MuMIn::r.squaredGLMM(m_base)[1, "R2m"]
r2_sal    <- MuMIn::r.squaredGLMM(m_sal)[1, "R2m"]
delta_r2  <- r2_sal - r2_base

lrt_sal   <- anova(m_base, m_sal)

effectsize_out <- data.frame(
  analysis         = "Salience on |dur_b1_vs_b4|",
  b_Salience       = sal_coef["Estimate"],
  se_Salience      = sal_coef["Std. Error"],
  t_Salience       = t_sal,
  df_resid         = df_resid,
  p_Salience       = sal_coef["Pr(>|t|)"],
  r_partial        = r_partial,
  eta2_partial     = eta2_part,
  R2m_base         = r2_base,
  R2m_with_sal     = r2_sal,
  delta_R2m        = delta_r2,
  lrt_chi2         = lrt_sal$Chisq[2],
  lrt_p            = lrt_sal$`Pr(>Chisq)`[2],
  n_trials         = nrow(qc),
  n_participants   = dplyr::n_distinct(qc$id)
)

# ── 3. Slope trajectory analysis ─────────────────────────────────────────────
# For each slope metric, fit the same structure as section 2 and extract
# the Salience coefficient. Tests whether Salience affects the *shape*
# of the breathing trajectory (not just total excursion).

slope_outcomes <- c("linear_slope", "early_change", "late_change", "mid_change")

slope_results <- purrr::map_dfr(slope_outcomes, function(out) {
  d <- qc |> dplyr::filter(!is.na(.data[[out]]))

  m_b <- lmerTest::lmer(
    reformulate(c("abs_delta", "(1 | id)"), response = out),
    data = d, REML = FALSE
  )
  m_s <- lmerTest::lmer(
    reformulate(c("salience", "abs_delta", "(1 | id)"), response = out),
    data = d, REML = FALSE
  )

  coef_row <- tryCatch(
    summary(m_s)$coefficients["salienceLow", ],
    error = function(e) NULL
  )
  if (is.null(coef_row)) return(NULL)

  t_v <- coef_row["t value"]
  df_v <- coef_row["df"]
  r_p  <- t_v / sqrt(t_v^2 + df_v)
  lrt  <- anova(m_b, m_s)

  data.frame(
    outcome     = out,
    b_Salience  = coef_row["Estimate"],
    se_Salience = coef_row["Std. Error"],
    t_Salience  = t_v,
    df_resid    = df_v,
    p_Salience  = coef_row["Pr(>|t|)"],
    r_partial   = r_p,
    lrt_chi2    = lrt$Chisq[2],
    lrt_p       = lrt$`Pr(>Chisq)`[2],
    n_trials    = nrow(d),
    n_ppts      = dplyr::n_distinct(d$id),
    row.names   = NULL
  )
})

cat("\n--- Slope trajectory: Salience effect ---\n")
print(slope_results, digits = 3, row.names = FALSE)


# ── 4. Salience × Accuracy on physio magnitude ───────────────────────────────
# Join behavioral Accuracy using content-based join:
#   id × Salience × Direction × row_number within cell
# (matches the strategy in utils.R::add_arousal_keyed)

behav_raw <- s5_long_breath

# Derive Salience and Direction in belt data from raw columns
qc_dir <- qc |>
  dplyr::mutate(
    Direction = dplyr::case_when(
      delta < 0 ~ "Faster",
      delta > 0 ~ "Slower",
      TRUE      ~ NA_character_
    ),
    Salience = dplyr::case_when(
      tolower(salience) == "high" ~ "High",
      tolower(salience) == "low"  ~ "Low",
      TRUE                        ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(Direction), !is.na(Salience))

qc_indexed <- qc_dir |>
  dplyr::group_by(id, Salience, Direction) |>
  dplyr::mutate(trial_in_cell = dplyr::row_number()) |>
  dplyr::ungroup()

behav_indexed <- behav_raw |>
  dplyr::filter(
    Condition == "breath",
    Direction %in% c("Faster", "Slower"),
    !is.na(Accuracy)
  ) |>
  dplyr::group_by(id, Salience, Direction) |>
  dplyr::mutate(trial_in_cell = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::select(id, Salience, Direction, trial_in_cell, Accuracy)

qc_acc <- qc_indexed |>
  dplyr::left_join(behav_indexed,
                   by = c("id", "Salience", "Direction", "trial_in_cell"))

match_pct <- mean(!is.na(qc_acc$Accuracy)) * 100
cat(sprintf("\nAccuracy join: %.1f%% trials matched (N=%d)\n",
            match_pct, sum(!is.na(qc_acc$Accuracy))))

# Stratified models: within Accuracy level, test Salience
acc_subsets <- list(
  correct   = dplyr::filter(qc_acc, Accuracy == 1),
  incorrect = dplyr::filter(qc_acc, Accuracy == 0),
  all       = qc_acc
)

acc_results <- purrr::map_dfr(names(acc_subsets), function(nm) {
  d <- acc_subsets[[nm]] |> dplyr::filter(!is.na(abs_dur_b1_vs_b4))
  if (nrow(d) < 50) return(NULL)

  m_b <- lmerTest::lmer(abs_dur_b1_vs_b4 ~ abs_delta + (1 | id),
                        data = d, REML = FALSE)
  m_s <- lmerTest::lmer(abs_dur_b1_vs_b4 ~ salience + abs_delta + (1 | id),
                        data = d, REML = FALSE)

  coef_row <- tryCatch(
    summary(m_s)$coefficients["salienceLow", ],
    error = function(e) NULL
  )
  if (is.null(coef_row)) return(NULL)

  t_v <- coef_row["t value"]
  df_v <- coef_row["df"]
  r_p  <- t_v / sqrt(t_v^2 + df_v)
  lrt  <- anova(m_b, m_s)

  data.frame(
    subset      = nm,
    b_Salience  = coef_row["Estimate"],
    se_Salience = coef_row["Std. Error"],
    t_Salience  = t_v,
    df_resid    = df_v,
    p_Salience  = coef_row["Pr(>|t|)"],
    r_partial   = r_p,
    lrt_chi2    = lrt$Chisq[2],
    lrt_p       = lrt$`Pr(>Chisq)`[2],
    n_trials    = nrow(d),
    n_ppts      = dplyr::n_distinct(d$id),
    row.names   = NULL
  )
})

cat("\n--- Salience × Accuracy: per-subset Salience effect ---\n")
print(acc_results, digits = 3, row.names = FALSE)

# Interaction model
qc_acc_matched <- qc_acc |> dplyr::filter(!is.na(Accuracy))
if (nrow(qc_acc_matched) > 100) {
  m_int <- lmerTest::lmer(
    abs_dur_b1_vs_b4 ~ salience * factor(Accuracy) + abs_delta + (1 | id),
    data = qc_acc_matched, REML = FALSE
  )
  cat("\n--- Salience × Accuracy interaction model ---\n")
  print(broom.mixed::tidy(m_int, effects = "fixed", conf.int = TRUE), n = Inf)
}


# Consolidated: table_belt_salience_followup.csv
# All three salience follow-up analyses in one file.
# component column identifies the analysis type.

table_belt_salience_followup <- dplyr::bind_rows(
  effectsize_out |>
    dplyr::mutate(component = "effect_size",
                  label     = analysis) |>
    dplyr::select(component, label, dplyr::everything(), -analysis),
  slope_results |>
    dplyr::mutate(component = "slope_trajectory",
                  label     = outcome) |>
    dplyr::rename(n_participants = n_ppts) |>
    dplyr::select(component, label, dplyr::everything(), -outcome),
  acc_results |>
    dplyr::mutate(component = "accuracy_stratified",
                  label     = subset) |>
    dplyr::rename(n_participants = n_ppts) |>
    dplyr::select(component, label, dplyr::everything(), -subset)
)

print(table_belt_salience_followup)
readr::write_csv(table_belt_salience_followup,
                 file.path(RESULTS_DIR, "table_belt_salience_followup.csv"))

cat("\nDone. Output written to:", RESULTS_DIR, "\n")
cat("  table_belt_salience_followup.csv\n")

