# ============================================================
# analysis_belt.R
# Study 5 belt physio analyses (primary + exploratory).
# Sources after analysis_arousal.R; uses s5l, s5s from MainAnalysis.R.
# Loads qcFile.xlsx (belt-specific; not in MainAnalysis.R).
#
# Sections:
#   1. Load data + regime subsets
#   2. Physio signal quality
#   3. Physio -> arousal (hit/miss; direction compliance)
#   F. Figures
#   S. Salience independence check
# Outputs: table_belt_physio_arousal.csv, table_belt_compliance.csv,
#          belt_direction_compliance_misses.csv, belt_salience_independence.csv,
#          belt_regime_key_terms.csv
# ============================================================

source(file.path(ANALYSIS_DIR, "theme_bcat.R"))

belt_fig_dir   <- file.path(FIG_DIR,   "Study5_Belt")
dir.create(belt_fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Belt-specific functions (moved from utils.R and models.R) ─

# ============================================================
# U4. Build four belt QC exclusion regime subsets from qcFull/qcSummary.
# `qcFull`:    FullResults sheet from qcFile.xlsx
# `qcSummary`: ResultsSummary sheet from qcFile.xlsx
#
# NB: belt_quality levels are "good", "degraded", "unusable" —
#   NOT "poor". Check your qcFile version if R3 is unexpectedly large.
# ============================================================
build_regime_datasets <- function(qcFull, qcSummary) {

  # ── Regime 2b: identify poor-synchrony participants ────────
  qcSummary <- qcSummary |>
    dplyr::mutate(
      breath_run = dplyr::if_else(first_condition == "Breath", 1L, 2L)
    )

  regime2b_exclude <- qcSummary |>
    dplyr::mutate(id = as.integer(id)) |>
    dplyr::filter(run == breath_run) |>
    dplyr::filter(
      r_delta_change_all  < 0.30,
      pct_correct_breaths < 50
    ) |>
    dplyr::pull(id)

  # ── Regime 3: unusable belt ─────────────────────────────────
  regime3_exclude <- qcSummary |>
    dplyr::mutate(id = as.integer(id)) |>
    dplyr::filter(belt_quality == "unusable") |>
    dplyr::pull(id) |>
    unique()

  message(sprintf(
    "Regime exclusions — R2b: %d participants  |  R3 additional: %d",
    length(regime2b_exclude),
    length(setdiff(regime3_exclude, regime2b_exclude))
  ))

  # ── Trial-level QC flags ────────────────────────────────────
  # The `condition` column in qcFull is contaminated with pipeline
  # logfile artefacts and cannot be used for filtering.
  # Reconstruct all classification variables from authoritative sources:
  #
  #   Condition:  run=1 & first_condition=="Visual" -> "visual", else "breath"
  #   Salience:   from `salience` column (capitalise first letter)
  #   Direction:  from sign of `delta`
  #               delta < 0 -> "Faster" (acceleration)
  #               delta > 0 -> "Slower" (deceleration)
  #               delta == 0 -> "NoChange"

  qc_all <- qcFull |>
    dplyr::filter(!is.na(id), !is.na(run)) |>
    dplyr::mutate(
      id = as.integer(id),
      # Reconstruct Condition
      Condition = dplyr::if_else(
        run == 1 & first_condition == "Visual",
        "visual", "breath"
      ),
      # Reconstruct Salience (standardise capitalisation)
      Salience = dplyr::case_when(
        tolower(salience) == "high" ~ "High",
        tolower(salience) == "low"  ~ "Low",
        TRUE                        ~ NA_character_
      ),
      # Reconstruct Direction from delta
      Direction = dplyr::case_when(
        delta < 0  ~ "Faster",
        delta > 0  ~ "Slower",
        delta == 0 ~ "NoChange",
        TRUE       ~ NA_character_
      )
    )

  # Paced change trials = breath condition, not NoChange
  qc_paced <- qc_all |>
    dplyr::filter(
      trial_available == TRUE,
      Condition       == "breath",
      Direction       != "NoChange",
      !is.na(Salience)
    ) |>
    dplyr::mutate(
      flag_cb     = correct_breaths == FALSE,
      flag_lag    = lag_flag        == TRUE,
      flag_dur    = is.na(dur_b1) | dur_b1 < 2.0 | dur_b1 > 7.0,
      n_flags     = as.integer(flag_cb) + as.integer(flag_lag) +
                    as.integer(flag_dur),
      regime2_bad = n_flags >= 2
    )

  # Sanity check: confirm Condition reconstruction
  cond_check <- qc_paced |>
    dplyr::count(first_condition, run, Condition)
  message("Condition reconstruction check (qc_paced):")
  print(as.data.frame(cond_check))

  # ── Four regime datasets ─────────────────────────────────────
  dat_r1  <- qc_paced
  dat_r2  <- dplyr::filter(qc_paced, !regime2_bad)
  dat_r2b <- dplyr::filter(qc_paced, !regime2_bad,
                            !id %in% regime2b_exclude)
  dat_r3  <- dplyr::filter(qc_paced, !regime2_bad,
                            !id %in% regime2b_exclude,
                            !id %in% regime3_exclude)

  regime_list <- list(R1 = dat_r1, R2 = dat_r2, R2b = dat_r2b, R3 = dat_r3)

  message("\nRegime sample sizes (paced trials):")
  for (rname in names(regime_list)) {
    d <- regime_list[[rname]]
    message(sprintf("  %s: %d trials, %d participants",
                    rname, nrow(d), dplyr::n_distinct(d$id)))
  }

  # Attach exclusion vectors as attributes for downstream use
  attr(regime_list, "regime2b_exclude") <- regime2b_exclude
  attr(regime_list, "regime3_exclude")  <- regime3_exclude

  regime_list
}


# ============================================================
# U5. add_arousal_keyed()
#
# Joins trial-level Arousal, Accuracy, Confidence, and Change
# from the staircase long file into a belt regime dataset.
#
# trial_num in qcFull is a global counter across all 80 trials
# (including NoChange). It is NOT sequential within change trials.
# The only reliable join strategy is content-based:
#   id × Salience × Direction × row_number within that cell
#
# Both qc_paced (belt) and long_paced (staircase) contain the
# same change trials in the same order within each cell.
#
# `dat`:        one regime dataset from build_regime_datasets()
# `long_paced`: staircase long data filtered to paced change trials,
#               with Salience and Direction columns present.
# ============================================================
add_arousal_keyed <- function(dat, long_paced) {

  # Index within id x Salience x Direction cell in long_paced
  long_indexed <- long_paced |>
    dplyr::group_by(id, Salience, Direction) |>
    dplyr::mutate(trial_in_cell = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::select(id, Salience, Direction, trial_in_cell,
                  Arousal, Accuracy, Confidence, Change)

  # Same index in the belt regime dataset
  dat_indexed <- dat |>
    dplyr::group_by(id, Salience, Direction) |>
    dplyr::mutate(trial_in_cell = dplyr::row_number()) |>
    dplyr::ungroup()

  dat_indexed |>
    dplyr::left_join(long_indexed,
                     by = c("id", "Salience", "Direction", "trial_in_cell"))
}

# ============================================================
# M6. Run H4 arousal models across all four belt exclusion regimes.
#   Names must be c("R1", "R2", "R2b", "R3").
# `group_var`:      passed to fit_arousal_models() for H4C.
#
# Returns:
#   $estimates  — tidy data frame of all fixed effects × regime
#   $key_terms  — subset to Change and Change:Accuracy rows
#   $models     — list of model-result lists per regime
# ============================================================
run_H4_regime_loop <- function(regime_arousal,
                               group_var   = NULL,
                               study_label = "Study5") {

  stopifnot(all(c("R1", "R2", "R2b", "R3") %in% names(regime_arousal)))

  all_estimates <- purrr::map_dfr(
    names(regime_arousal),
    function(rname) {
      dat <- regime_arousal[[rname]] |>
        dplyr::filter(!is.na(Arousal), !is.na(Change)) |>
        dplyr::mutate(Change2 = Change^2)   # pre-compute — Change^2 in formula is wrong

      cat(sprintf("\n  Running regime %s: %d trials, %d participants\n",
                  rname, nrow(dat), dplyr::n_distinct(dat$id)))

      # Fit H4A and H4B
      mA <- lmerTest::lmer(Arousal ~ Change + Change2 + (1 | id),
                        data = dat, REML = FALSE)
      mB <- lmerTest::lmer(
        Arousal ~ Change * Accuracy + Change2 + (1 | id),
        data = dat |> dplyr::filter(!is.na(Accuracy)),
        REML = FALSE
      )

      dplyr::bind_rows(
        broom.mixed::tidy(mA, effects = "fixed", conf.int = TRUE) |>
          dplyr::mutate(hypothesis = "H4A", regime = rname),
        broom.mixed::tidy(mB, effects = "fixed", conf.int = TRUE) |>
          dplyr::mutate(hypothesis = "H4B", regime = rname)
      )
    }
  )

  key_terms <- all_estimates |>
    dplyr::filter(
      (hypothesis == "H4A" & term == "Change") |
      (hypothesis == "H4B" & term %in% c("Change:Accuracy", "Accuracy:Change"))
    ) |>
    dplyr::mutate(
      regime = factor(regime, levels = c("R1", "R2", "R2b", "R3")),
      term_clean = dplyr::if_else(
        term %in% c("Change:Accuracy", "Accuracy:Change"),
        "Change:Accuracy", term
      ),
      term_label = dplyr::recode(
        term_clean,
        "Change"          = "H4A: Change -> Arousal",
        "Change:Accuracy" = "H4B: Change x Accuracy"
      )
    )

  cat(sprintf(
    "\n[%s] Regime comparison — H4B Change×Accuracy across regimes:\n",
    study_label
  ))
  print(
    key_terms |>
      dplyr::filter(term == "Change:Accuracy") |>
      dplyr::select(regime, estimate, std.error, statistic,
                    p.value, conf.low, conf.high) |>
      as.data.frame(),
    digits = 3, row.names = FALSE
  )

  list(estimates = all_estimates, key_terms = key_terms)
}

# ============================================================
# 1. LOAD DATA
# ============================================================
message("Loading belt-specific data...")

# s5l, s5s, qcFull, qcSummary all loaded by MainAnalysis.R
long         <- s5l
summary_data <- s5s

# ── Prepare paced staircase long for arousal join ────────────
# Includes both sessions; all breath-condition paced trials.
paced_conds <- c("highSalienceAcc", "highSalienceDec",
                 "lowSalienceAcc",  "lowSalienceDec")

# Restrict to breath-condition trials only and sort by trial onset time
# before indexing. Visual-condition trials must be excluded: including them
# creates index collisions within id x Salience x Direction that corrupt
# the Accuracy join (visual and breath trials get the same trial_in_cell).
long_paced <- long |>
  dplyr::filter(
    taskCondition %in% paced_conds,
    Condition == "breath"
  ) |>
  dplyr::mutate(id = as.integer(id)) |>
  dplyr::arrange(id, Salience, Direction, trial.started) |>
  dplyr::group_by(id, Salience, Direction) |>
  dplyr::mutate(trial_in_cell = dplyr::row_number()) |>
  dplyr::ungroup()

message(sprintf("Paced staircase trials (both sessions): %d across %d participants",
                nrow(long_paced), dplyr::n_distinct(long_paced$id)))


# ============================================================
# 2. BUILD REGIME DATASETS
# ============================================================
regimes       <- build_regime_datasets(qcFull, qcSummary)
regime_arousal <- lapply(regimes, add_arousal_keyed, long_paced = long_paced)

# Join quality summary for exploratory belt models
qcSummary_breath <- qcSummary |>
  dplyr::mutate(
    breath_run = dplyr::if_else(first_condition == "Breath", 1L, 2L)
  ) |>
  dplyr::filter(run == breath_run) |>
  dplyr::select(id, pct_correct_breaths, r_delta_change_all,
                belt_quality, median_correlation)

# Spot-check join quality
n_joined <- sum(!is.na(regime_arousal$R1$Arousal))
n_total  <- nrow(regime_arousal$R1)
message(sprintf("Arousal join (R1): %d / %d trials matched (%.0f%%)",
                n_joined, n_total, 100 * n_joined / n_total))


# ============================================================
# 3. SIGN ALIGNMENT DIAGNOSTIC
#
# Verify dur_b1_vs_b4 and Change share the same sign convention
# before any physio models are run. Positive dur_b1_vs_b4 means
# last breath longer than first (breathed slower), corresponding
# to positive Change. Expect a positive correlation on paced trials.
# ============================================================
message("\n--- Sign alignment diagnostic ---")

sign_check <- regime_arousal$R1 |>
  dplyr::filter(!is.na(dur_b1_vs_b4), !is.na(delta)) |>
  dplyr::summarise(
    r_overall  = cor(dur_b1_vs_b4, delta, use = "complete.obs"),
    r_faster   = cor(dur_b1_vs_b4[delta < 0],
                     delta[delta < 0], use = "complete.obs"),
    r_slower   = cor(dur_b1_vs_b4[delta > 0],
                     delta[delta > 0], use = "complete.obs"),
    n          = sum(!is.na(dur_b1_vs_b4) & !is.na(delta))
  )

message(sprintf(
  "  r(dur_b1_vs_b4, intended_change) overall = %.3f  faster = %.3f  slower = %.3f",
  sign_check$r_overall, sign_check$r_faster, sign_check$r_slower
))

if (sign_check$r_overall < 0) {
  stop("SIGN MISMATCH: dur_b1_vs_b4 is negatively correlated with intended ",
       "change. Check qcFull variable definition before proceeding.")
} else {
  message("  Sign alignment confirmed — positive correlation as expected.")
}

# Use Change (canonical, from long_paced join) as intended change predictor.
# dur_b1_vs_b4 from qcFull is the observed physio predictor.


# ============================================================
# LAYER 1: REGIME SENSITIVITY ANALYSIS (preregistered)
# ============================================================
cat("\n\n========================================\n")
cat("LAYER 1: REGIME SENSITIVITY (H4A, H4B)\n")
cat("Preregistered sensitivity analysis\n")
cat("========================================\n")

# H4C is behavioural only — visual-condition trials are not in qc_paced.
# Do NOT include visual-condition trials in regime analyses.
# H4C is run separately on long in unified_analysis.R.

regime_results <- run_H4_regime_loop(
  regime_arousal = regime_arousal,
  study_label    = "Study5"
)

# Export regime key terms for standalone figure script and reproducibility
readr::write_csv(
  regime_results$key_terms,
  file.path(RESULTS_DIR, "belt_regime_key_terms.csv")
)
message("Saved: belt_regime_key_terms.csv")

# ── Extract primary (R1) p-values for BH correction ─────────
h4a_r1_p <- regime_results$estimates |>
  dplyr::filter(regime == "R1", hypothesis == "H4A", term == "Change") |>
  dplyr::pull(p.value)

h4b_r1_p <- regime_results$estimates |>
  dplyr::filter(regime == "R1", hypothesis == "H4B",
                stringr::str_detect(term, ":")) |>
  dplyr::pull(p.value)

message(sprintf("\nR1 primary results: H4A p=%.4f  H4B p=%.4f",
                h4a_r1_p, h4b_r1_p))

# ── Regime agreement summary ─────────────────────────────────
regime_agreement <- regime_results$key_terms |>
  dplyr::group_by(term_label) |>
  dplyr::summarise(
    direction_consistent = all(estimate > 0) || all(estimate < 0),
    min_estimate = min(estimate),
    max_estimate = max(estimate),
    range        = max(estimate) - min(estimate),
    all_sig      = all(p.value < .05),
    .groups = "drop"
  )

cat("\nRegime consistency summary:\n")
print(as.data.frame(regime_agreement), row.names = FALSE)


# ============================================================
# LAYER 2: PHYSIOLOGY AS PREDICTOR (exploratory)
# ============================================================
cat("\n\n========================================\n")
cat("LAYER 2: PHYSIOLOGY AS AROUSAL PREDICTOR\n")
cat("(Exploratory — not preregistered)\n")
cat("========================================\n")

# Use Regime 2 throughout LAYER 2 and 3:
# trial-level QC applied (n_flags < 2), giving cleaner physio signal
# while retaining reasonable sample size.
explore_base <- regime_arousal$R2 |>
  dplyr::filter(!is.na(Arousal), !is.na(dur_b1_vs_b4), !is.na(Change)) |>
  dplyr::mutate(Change2 = Change^2)


# ── 2a. Simple physio → arousal link ─────────────────────────
cat("\n--- 2a. Arousal ~ dur_b1_vs_b4 (physio only) ---\n")
m2a <- lmerTest::lmer(Arousal ~ dur_b1_vs_b4 + (1 | id),
                   data = explore_base, REML = FALSE)
print(broom.mixed::tidy(m2a, effects = "fixed", conf.int = TRUE), n = Inf)
cat(sprintf("  R²m=%.3f  R²c=%.3f\n",
            MuMIn::r.squaredGLMM(m2a)[1], MuMIn::r.squaredGLMM(m2a)[2]))


# ── 2b. Physio beyond intended change ────────────────────────
# Key test: does observed compliance add unique predictive value
# over and above the experimental manipulation?
cat("\n--- 2b. Arousal ~ Change + dur_b1_vs_b4 + Change² ---\n")

m2b_base <- lmerTest::lmer(Arousal ~ Change + Change2 + (1 | id),
                        data = explore_base, REML = FALSE)
m2b_full <- lmerTest::lmer(Arousal ~ Change + Change2 + dur_b1_vs_b4 + (1 | id),
                        data = explore_base, REML = FALSE)

lrt_2b <- anova(m2b_base, m2b_full)
cat("  LRT (does dur_b1_vs_b4 add beyond Change?):\n")
print(lrt_2b)
print(broom.mixed::tidy(m2b_full, effects = "fixed", conf.int = TRUE), n = Inf)
cat(sprintf("  R²m base=%.3f  full=%.3f  dR2m=%.3f\n",
            MuMIn::r.squaredGLMM(m2b_base)[1],
            MuMIn::r.squaredGLMM(m2b_full)[1],
            MuMIn::r.squaredGLMM(m2b_full)[1] -
              MuMIn::r.squaredGLMM(m2b_base)[1]))


# ── 2c. Belt quality as moderator ────────────────────────────
# Does the strength of the physio → arousal link scale with
# how well participants complied (continuous measure)?
cat("\n--- 2c. Belt quality as continuous moderator ---\n")

explore_qual <- explore_base |>
  dplyr::left_join(qcSummary_breath |>
                     dplyr::select(id, pct_correct_breaths,
                                   r_delta_change_all),
                   by = "id") |>
  dplyr::filter(!is.na(pct_correct_breaths)) |>
  dplyr::mutate(
    pct_correct_z = scale(pct_correct_breaths)[, 1]
  )

m2c <- lmerTest::lmer(
  Arousal ~ Change * pct_correct_z + Change2 + (1 | id),
  data = explore_qual, REML = FALSE
)
print(broom.mixed::tidy(m2c, effects = "fixed", conf.int = TRUE), n = Inf)


# ============================================================
# LAYER 2d: DECOMPOSITION MODEL
#
# Partitions arousal variance into three orthogonal components:
#   (1) Intended change  — the experimental manipulation
#   (2) Observed compliance — physio beyond what was demanded
#       (residualised on Change to make it orthogonal)
#   (3) Conscious detection — Accuracy
#
# Nested LRT sequence reveals unique contribution of each layer.
# Standardised predictors allow direct β comparison.
# ============================================================
cat("\n--- 2d. Decomposition: intended change vs compliance vs detection ---\n")

explore_decomp <- explore_base |>
  dplyr::filter(!is.na(Accuracy)) |>
  dplyr::mutate(
    Change_z        = scale(Change)[, 1],
    # Residualise dur_b1_vs_b4 on Change so the two physio predictors
    # are orthogonal — separates what the instruction demanded from
    # how the participant's body actually responded beyond that
    dur_resid       = residuals(lm(dur_b1_vs_b4 ~ Change,
                                   data = dplyr::pick(everything()))),
    dur_resid_z     = scale(dur_resid)[, 1]
  )

cat(sprintf("  Decomposition N: %d trials, %d participants\n",
            nrow(explore_decomp), dplyr::n_distinct(explore_decomp$id)))
cat(sprintf("  r(Change, dur_resid) after orthogonalisation: %.4f (expect ~0)\n",
            cor(explore_decomp$Change_z, explore_decomp$dur_resid_z)))

# M1: intended change only
m_d1 <- lmerTest::lmer(Arousal ~ Change_z + (1 | id),
                    data = explore_decomp, REML = FALSE)

# M2: + observed compliance (orthogonal to intended change)
m_d2 <- lmerTest::lmer(Arousal ~ Change_z + dur_resid_z + (1 | id),
                    data = explore_decomp, REML = FALSE)

# M3: + conscious detection
m_d3 <- lmerTest::lmer(Arousal ~ Change_z + dur_resid_z + Accuracy + (1 | id),
                    data = explore_decomp, REML = FALSE)

# M4: + compliance × detection interaction
# Key test: does detection amplify the physio → arousal link,
# or does the physio pathway operate independently of awareness?
m_d4 <- lmerTest::lmer(
  Arousal ~ Change_z + dur_resid_z * Accuracy + (1 | id),
  data = explore_decomp, REML = FALSE
)

decomp_models <- list(
  "M1: Intended change"              = m_d1,
  "M2: + Observed compliance"        = m_d2,
  "M3: + Conscious detection"        = m_d3,
  "M4: + Compliance × Detection"     = m_d4
)

cat("\n  Model comparison:\n")
for (nm in names(decomp_models)) {
  m   <- decomp_models[[nm]]
  r2  <- MuMIn::r.squaredGLMM(m)
  cat(sprintf("  %-35s  AIC=%7.1f  R²m=%.3f  R²c=%.3f\n",
              nm, AIC(m), r2[1, "R2m"], r2[1, "R2c"]))
}

cat("\n  LRT sequence M1→M2→M3→M4:\n")
print(anova(m_d1, m_d2, m_d3, m_d4))

cat("\n  Full model (M4) fixed effects:\n")
print(broom.mixed::tidy(m_d4, effects = "fixed", conf.int = TRUE), n = Inf)

cat("\n  Interpretation:\n")
cat("  M1→M2: does observed compliance add beyond intended change?\n")
cat("  M2→M3: does conscious detection add beyond physio compliance?\n")
cat("  M3→M4: does detection amplify the compliance→arousal link?\n")
cat("  Null M3→M4 interaction = physio pathway independent of awareness\n")


# ============================================================
# LAYER 3: NON-CONSCIOUS COMPLIANCE (named exploratory centrepiece)
#
# On missed trials (Accuracy = 0):
#   3a. Binomial test — is direction_correct above chance (50%)?
#   3b. Mixed model   — is belt change magnitude similar to hits?
#   3c. Multilevel    — does physio predict arousal on missed trials?
#        Person-mean centering decomposes within- vs between-person
#        physio-arousal covariation in a single multilevel model.
#   3d. Interaction   — does detection moderate the physio→arousal link?
#        Theory predicts NO — misattribution operates without awareness.
# ============================================================
cat("\n\n========================================\n")
cat("LAYER 3: NON-CONSCIOUS COMPLIANCE\n")
cat("(Named exploratory — not preregistered)\n")
cat("========================================\n")

# Change trials only (exclude NoChange — no direction signal)
explore_change <- regime_arousal$R2 |>
  dplyr::filter(
    !is.na(Accuracy),
    !is.na(direction_correct),
    !is.na(dur_b1_vs_b4),
    delta != 0   # change trials only
  ) |>
  dplyr::mutate(Change2 = Change^2)

explore_miss <- dplyr::filter(explore_change, Accuracy == 0)
explore_hit  <- dplyr::filter(explore_change, Accuracy == 1)

message(sprintf("Change trials (R2): %d total  |  hits=%d  misses=%d",
                nrow(explore_change),
                nrow(explore_hit),
                nrow(explore_miss)))


# ── 3a. Direction compliance on missed trials ─────────────────
cat("\n--- 3a. Direction compliance by detection outcome ---\n")

dir_summary <- explore_change |>
  dplyr::group_by(Accuracy) |>
  dplyr::summarise(
    n               = dplyr::n(),
    pct_dir_correct = mean(direction_correct, na.rm = TRUE) * 100,
    mean_belt_mag   = mean(abs(dur_b1_vs_b4), na.rm = TRUE),
    sd_belt_mag     = sd(abs(dur_b1_vs_b4),   na.rm = TRUE),
    .groups = "drop"
  )

cat("Direction correct rate by detection accuracy:\n")
print(as.data.frame(dir_summary), digits = 3, row.names = FALSE)

# Binomial test: is direction compliance on misses above chance?
n_miss_dir  <- sum(!is.na(explore_miss$direction_correct))
k_miss_dir  <- sum(explore_miss$direction_correct == 1, na.rm = TRUE)
binom_miss  <- binom.test(k_miss_dir, n_miss_dir, p = 0.5)

cat(sprintf("\n  Binomial test (misses vs 50%% chance):\n"))
cat(sprintf("  Direction correct on misses: %d / %d = %.1f%%\n",
            k_miss_dir, n_miss_dir,
            100 * k_miss_dir / n_miss_dir))
cat(sprintf("  p = %.4f  95%% CI [%.3f, %.3f]\n",
            binom_miss$p.value,
            binom_miss$conf.int[1],
            binom_miss$conf.int[2]))


# ── 3b. Belt magnitude hits vs misses ─────────────────────────
# If similar magnitude on misses and hits, non-detection is a
# failure of conscious registration, not of physical compliance.
cat("\n--- 3b. Belt change magnitude: hits vs misses ---\n")

m3b <- lmerTest::lmer(
  abs(dur_b1_vs_b4) ~ Accuracy + (1 | id),
  data = explore_change, REML = FALSE
)
cat(sprintf("  Accuracy effect on |dur_b1_vs_b4|: %s\n",
            fmt_lmer(m3b, "Accuracy")))

# Non-parametric check (Wilcoxon on trial-level — not nested)
wt_belt <- wilcox.test(
  explore_hit$dur_b1_vs_b4,
  explore_miss$dur_b1_vs_b4
)
cat(sprintf("  Wilcoxon W = %.0f  p = %.4f\n",
            wt_belt$statistic, wt_belt$p.value))
cat("  (Non-nested Wilcoxon for descriptive check only; LMM is primary)\n")

# Save 3a and 3b results to CSV
m3b_cf <- summary(m3b)$coefficients
table_belt_compliance <- tibble::tibble(
  # 3a: direction compliance on missed trials
  n_miss_dir          = n_miss_dir,
  k_miss_dir_correct  = k_miss_dir,
  pct_miss_dir_correct = round(100 * k_miss_dir / n_miss_dir, 2),
  pct_hit_dir_correct  = round(dir_summary$pct_dir_correct[dir_summary$Accuracy == 1], 2),
  binom_p             = binom_miss$p.value,
  binom_ci_lower      = binom_miss$conf.int[1],
  binom_ci_upper      = binom_miss$conf.int[2],
  # 3b: belt magnitude hits vs misses
  b_Accuracy_3b       = m3b_cf["Accuracy", "Estimate"],
  se_Accuracy_3b      = m3b_cf["Accuracy", "Std. Error"],
  t_Accuracy_3b       = m3b_cf["Accuracy", "t value"],
  p_Accuracy_3b       = m3b_cf["Accuracy", "Pr(>|t|)"],
  wilcox_W            = unname(wt_belt$statistic),
  wilcox_p            = wt_belt$p.value
)
readr::write_csv(table_belt_compliance,
                 file.path(RESULTS_DIR, "table_belt_compliance.csv"))
message("Saved: table_belt_compliance.csv")


# ── 3c. Physio → Arousal on missed trials: multilevel model ───
#
# Person-mean centering decomposes dur_b1_vs_b4 into:
#   WITHIN:   trial-level deviations from participant's own mean
#             (does THIS trial's physio predict THIS trial's arousal?)
#   BETWEEN:  participant's mean physio across all missed trials
#             (do more compliant participants report higher arousal?)
#
# Theory predicts the WITHIN-person effect should be significant —
# misattribution is a dynamic trial-by-trial mechanism.
# If only BETWEEN is significant, the link is a stable trait, not
# a misattribution process.
#
# We also compute the person-mean across ALL trials (not just misses)
# as a Level-2 covariate, to prevent the between-miss term from
# soaking up general arousal-compliance trait variance.
# ============================================================
cat("\n--- 3c. Person-mean-centred physio → arousal on missed trials ---\n")

# Person-means across ALL change trials (Level-2 covariate)
person_means_all <- explore_change |>
  dplyr::group_by(id) |>
  dplyr::summarise(
    dur_mean_all = mean(dur_b1_vs_b4, na.rm = TRUE),
    .groups = "drop"
  )

# Person-means across MISSED trials only
person_means_miss <- explore_miss |>
  dplyr::group_by(id) |>
  dplyr::summarise(
    dur_mean_miss = mean(dur_b1_vs_b4, na.rm = TRUE),
    .groups = "drop"
  )

# Prepare missed-trial data with centred predictors
miss_centred <- explore_miss |>
  dplyr::left_join(person_means_miss, by = "id") |>
  dplyr::left_join(person_means_all,  by = "id") |>
  dplyr::filter(!is.na(Arousal)) |>
  dplyr::mutate(
    # Within: trial deviation from person's own miss-trial mean
    dur_within  = dur_b1_vs_b4 - dur_mean_miss,
    # Between (misses): person's mean on missed trials
    dur_between = dur_mean_miss,
    # Level-2 control: person's mean across ALL trials
    dur_all_c   = scale(dur_mean_all)[, 1]
  )

cat(sprintf("  Missed-trial N: %d trials, %d participants\n",
            nrow(miss_centred), dplyr::n_distinct(miss_centred$id)))

# Primary multilevel model: within + between decomposition
m3c_main <- lmerTest::lmer(
  Arousal ~ dur_within + dur_between + Change2 + (1 | id),
  data = miss_centred, REML = FALSE
)

# Robustness: control for compliance trait across all trials
m3c_robust <- lmerTest::lmer(
  Arousal ~ dur_within + dur_between + dur_all_c + Change2 + (1 | id),
  data = miss_centred, REML = FALSE
)

cat("\n  Primary model (within + between):\n")
print(broom.mixed::tidy(m3c_main, effects = "fixed", conf.int = TRUE),
      n = Inf)
cat(sprintf("  R²m=%.3f  R²c=%.3f\n",
            MuMIn::r.squaredGLMM(m3c_main)[1],
            MuMIn::r.squaredGLMM(m3c_main)[2]))

cat("\n  Robustness (controlling for overall compliance trait):\n")
print(broom.mixed::tidy(m3c_robust, effects = "fixed", conf.int = TRUE),
      n = Inf)

# Compare within-person slope on missed vs hit trials
# (same model structure on hits, for contrast)
hit_centred <- explore_hit |>
  dplyr::group_by(id) |>
  dplyr::mutate(dur_within = dur_b1_vs_b4 - mean(dur_b1_vs_b4, na.rm = TRUE),
                dur_between = mean(dur_b1_vs_b4, na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::filter(!is.na(Arousal)) |>
  dplyr::mutate(Change2 = Change^2)

m3c_hits <- lmerTest::lmer(
  Arousal ~ dur_within + dur_between + Change2 + (1 | id),
  data = hit_centred, REML = FALSE
)

# Extract within-person slopes for comparison
# Note: check actual term name from tidy output in case lmerTest renames it
.tidy_main <- broom.mixed::tidy(m3c_main, effects = "fixed", conf.int = TRUE)
.tidy_hits <- broom.mixed::tidy(m3c_hits, effects = "fixed", conf.int = TRUE)
message("  m3c_main terms: ", paste(.tidy_main$term, collapse = ", "))
message("  m3c_hits terms: ", paste(.tidy_hits$term, collapse = ", "))

# Use grep to find the within term regardless of exact name
.within_term <- grep("within|dur_within", .tidy_main$term, value = TRUE)[1]
message("  Matched within term: ", .within_term)

within_miss <- .tidy_main |>
  dplyr::filter(term == .within_term) |>
  dplyr::mutate(condition = "Misses")

within_hits <- .tidy_hits |>
  dplyr::filter(term == .within_term) |>
  dplyr::mutate(condition = "Hits")

within_comparison <- dplyr::bind_rows(within_hits, within_miss)

cat("\n  Within-person physio → arousal slope comparison:\n")
print(
  within_comparison |>
    dplyr::select(dplyr::any_of(
      c("condition", "estimate", "std.error", "statistic",
        "p.value", "conf.low", "conf.high")
    )) |>
    as.data.frame(),
  digits = 3, row.names = FALSE
)
cat("  (Similar slopes = physio→arousal pathway doesn't require awareness)\n")


# ── 3d. Formal interaction test ────────────────────────────────
# `Arousal ~ dur_b1_vs_b4 * Accuracy + Change² + (1|id)` on all change trials.
# Interaction = does detection amplify the physio→arousal link?
# Theory predicts a NULL interaction (misattribution independent of awareness).
# Support for null via Bayesian approach: run with BayesFactor or report
# the interaction LRT alongside the within-person model above.
cat("\n--- 3d. Formal interaction test: physio × Accuracy ---\n")

# Use within-person centred physio for the interaction model
explore_change_c <- explore_change |>
  dplyr::group_by(id) |>
  dplyr::mutate(
    dur_within  = dur_b1_vs_b4 - mean(dur_b1_vs_b4, na.rm = TRUE),
    dur_between = mean(dur_b1_vs_b4, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(!is.na(Arousal))

m3d_base <- lmerTest::lmer(
  Arousal ~ dur_within + dur_between + Accuracy + Change2 + (1 | id),
  data = explore_change_c, REML = FALSE
)

m3d_int  <- lmerTest::lmer(
  Arousal ~ dur_within * Accuracy + dur_between + Change2 + (1 | id),
  data = explore_change_c, REML = FALSE
)

lrt_3d <- anova(m3d_base, m3d_int)

cat("  LRT: does Accuracy moderate the within-person physio → arousal link?\n")
print(lrt_3d)
cat("\n  Full interaction model fixed effects:\n")
print(broom.mixed::tidy(m3d_int, effects = "fixed", conf.int = TRUE),
      n = Inf)
cat("\n  Interpretation key:\n")
cat("  dur_within main effect = physio→arousal WITHOUT detection\n")
cat("  dur_within:Accuracy    = ADDITIONAL amplification WHEN detected\n")
cat("  Null/small interaction → physio pathway operates without awareness\n")

# Export: table_belt_physio_arousal.csv
# Layer 3c: within-person physio -> arousal on missed vs hit trials.
# Layer 3d: formal physio x Accuracy interaction LRT.
# These are main-text results (belt corroboration of awareness gating).

.tidy_3d <- broom.mixed::tidy(m3d_int, effects = "fixed", conf.int = TRUE)
.int_term_3d <- grep("dur_within.*Accuracy|Accuracy.*dur_within",
                     .tidy_3d$term, value = TRUE)[1]

table_belt_physio_arousal <- dplyr::bind_rows(
  # 3c: within-person slope on missed vs hit trials
  within_comparison |>
    dplyr::select(dplyr::any_of(
      c("condition", "estimate", "std.error", "statistic",
        "p.value", "conf.low", "conf.high")
    )) |>
    dplyr::mutate(
      layer = "3c_physio_arousal_by_accuracy",
      note  = "Within-person dur_b1_vs_b4 -> Arousal; Misses = primary null test"),
  # 3d: interaction coefficient and LRT
  tibble::tibble(
    condition  = "dur_within x Accuracy (interaction)",
    estimate   = if (!is.na(.int_term_3d))
                   .tidy_3d$estimate[.tidy_3d$term == .int_term_3d]
                 else NA_real_,
    std.error  = if (!is.na(.int_term_3d))
                   .tidy_3d$std.error[.tidy_3d$term == .int_term_3d]
                 else NA_real_,
    p.value    = if (!is.na(.int_term_3d))
                   .tidy_3d$p.value[.tidy_3d$term == .int_term_3d]
                 else NA_real_,
    lrt_chi2   = lrt_3d$Chisq[2],
    lrt_p      = lrt_3d$`Pr(>Chisq)`[2],
    layer      = "3d_physio_x_accuracy_interaction",
    note       = "LRT: does detection moderate physio->arousal? Positive = amplification on hits"
  )
)

readr::write_csv(table_belt_physio_arousal,
                 file.path(RESULTS_DIR, "table_belt_physio_arousal.csv"))
message("Saved: table_belt_physio_arousal.csv")


# ── 3e. Direction compliance × arousal on missed trials ───────
#
# Within missed trials: does physically complying with the direction
# instruction (direction_correct == TRUE) predict higher arousal
# than non-compliant missed trials (direction_correct == FALSE)?
#
# This tests the minimum interoceptive registration hypothesis:
# if the body moves in the right direction but the participant
# doesn't consciously detect the change, does that directional
# compliance carry any arousal signal?
#
# If direction_correct predicts arousal on misses:
#   → the body's directional compliance is a carrier of the
#     nonconscious effect, not just any physiological change
# If not:
#   → nonconscious transfer (if present in 3c) is driven by
#     magnitude (dur_b1_vs_b4), not directional compliance per se
cat("\n--- 3e. Direction compliance × arousal on missed trials ---\n")

explore_miss_dir <- explore_miss |>
  dplyr::filter(!is.na(direction_correct), !is.na(Arousal)) |>
  dplyr::mutate(
    dir_correct_f = factor(direction_correct,
                            levels = c(0, 1),
                            labels = c("Non-compliant", "Compliant"))
  )

# Descriptive
dir_miss_desc <- explore_miss_dir |>
  dplyr::group_by(dir_correct_f) |>
  dplyr::summarise(
    M_arousal = mean(Arousal, na.rm = TRUE),
    SD_arousal = sd(Arousal,  na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )
cat("Arousal by direction compliance on missed trials:\n")
print(as.data.frame(dir_miss_desc), digits = 3, row.names = FALSE)

# Mixed model: direction_correct predicts arousal on misses?
m3e <- lmerTest::lmer(
  Arousal ~ direction_correct + abs(dur_b1_vs_b4) + Change2 + (1 | id),
  data = explore_miss_dir, REML = FALSE
)
cat("\nModel: Arousal ~ direction_correct + |physio| + Change² + (1|id)\n")
print(broom.mixed::tidy(m3e, effects = "fixed", conf.int = TRUE), n = Inf)

# Bayesian null for direction_correct effect
dir_miss_bf <- tryCatch(
  BayesFactor::ttestBF(
    x = explore_miss_dir$Arousal[explore_miss_dir$direction_correct == 1],
    y = explore_miss_dir$Arousal[explore_miss_dir$direction_correct == 0]
  ),
  error = function(e) NULL
)
if (!is.null(dir_miss_bf)) {
  bf10 <- BayesFactor::extractBF(dir_miss_bf)$bf
  cat(sprintf("\nBF10 = %.3f  BF01 = %.3f — direction compliance × arousal on misses\n",
              bf10, 1/bf10))
}

# Save
# direction_correct is logical in qcFile, so R codes the term as
# "direction_correctTRUE" (dummy-coded TRUE vs FALSE reference).
tidy_m3e <- broom.mixed::tidy(m3e, effects = "fixed")
dc_term   <- grep("direction_correct", tidy_m3e$term, value = TRUE)[1]
dir_miss_results <- tibble::tibble(
  analysis      = "direction_correct on missed trial arousal",
  b_dir_correct = tidy_m3e |>
    dplyr::filter(term == dc_term) |> dplyr::pull(estimate),
  p_dir_correct = tidy_m3e |>
    dplyr::filter(term == dc_term) |> dplyr::pull(p.value),
  BF10          = if (!is.null(dir_miss_bf))
    unname(BayesFactor::extractBF(dir_miss_bf)$bf) else NA_real_
)
readr::write_csv(dir_miss_results,
                 file.path(RESULTS_DIR, "belt_direction_compliance_misses.csv"))
message("Saved: belt_direction_compliance_misses.csv")
# ============================================================
cat("\n\n========================================\n")
cat("FIGURES\n")
cat("========================================\n")

# ── F1: Regime comparison forest plot (restyled) ──────────────
#
# Belt QC exclusion regimes:
#   R1  = all paced trials (primary, no exclusions)
#   R2  = trial-level QC applied (n_flags < 2)
#   R2b = R2 + poor-synchrony participants excluded
#   R3  = R2b + unusable-belt participants excluded
#
# Shared x-axis across H4A and H4B panels for direct comparison.
# Regime colours: dark-to-light gradient (R1 darkest = primary).
if (nrow(regime_results$key_terms) == 0) {
  message("  [WARNING] key_terms is empty — check that 'Change' and 'Change:Accuracy'",
          " terms were found in regime models. Skipping F1.")
} else {

  regime_colours <- c(
    R1  = "#2d3748",   # darkest  — primary result
    R2  = "#4a5568",
    R2b = "#718096",
    R3  = "#a0aec0"    # lightest — most restrictive
  )

  f1 <- regime_results$key_terms |>
    dplyr::mutate(
      regime = factor(regime, levels = c("R3", "R2b", "R2", "R1")),
      pt_size = dplyr::if_else(regime == "R1", 4, 3)
    ) |>
    ggplot2::ggplot(ggplot2::aes(
      x      = estimate,
      xmin   = conf.low,
      xmax   = conf.high,
      y      = regime,
      colour = regime
    )) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype   = "dashed",
      colour     = "grey55",
      linewidth  = 0.5
    ) +
    ggplot2::geom_errorbarh(height = 0.22, linewidth = 0.8) +
    ggplot2::geom_point(
      ggplot2::aes(size = regime == "R1")
    ) +
    ggplot2::scale_colour_manual(
      values = regime_colours,
      guide  = "none"
    ) +
    ggplot2::scale_size_manual(
      values = c("TRUE" = 4, "FALSE" = 3),
      guide  = "none"
    ) +
    # Fixed (shared) x-axis: H4A and H4B on identical scale
    ggplot2::facet_wrap(~ term_label, scales = "fixed", nrow = 1) +
    ggplot2::labs(
      title    = "H4 stability across belt QC exclusion regimes",
      subtitle = paste0(
        "\u03b2 \u00b1 95% CI  |  ",
        "R1 = primary (no QC exclusions);  ",
        "R2 = trial-level QC;  ",
        "R2b = + poor synchrony excluded;  ",
        "R3 = + unusable belt excluded"
      ),
      x = "\u03b2 (effect on Arousal)",
      y = "Exclusion regime"
    ) +
    theme_bcat(base_size = 12) +
    ggplot2::theme(
      plot.subtitle      = ggplot2::element_text(size = 8.5, colour = "grey45",
                                                  lineheight = 1.3),
      strip.text         = ggplot2::element_text(face = "bold", size = 11),
      axis.title.y       = ggplot2::element_blank(),
      axis.text.y        = ggplot2::element_text(size = 10),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "grey88",
                                                  linewidth = 0.3),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.spacing      = ggplot2::unit(2, "lines")
    )

  save_bcat_fig(file.path(belt_fig_dir, "S5_regime_comparison"),
                f1, width = 10, height = 3.8)
  ggplot2::ggsave(
    file.path(belt_fig_dir, "S5_regime_comparison.png"),
    f1, width = 10, height = 3.8, dpi = 300
  )
  message("Saved: S5_regime_comparison.png")
}


# ── F2: Decomposition model Decomposition model dR2m bar chart ────────────────────
decomp_r2 <- tibble::tibble(
  model = names(decomp_models),
  R2m   = sapply(decomp_models,
                 function(m) MuMIn::r.squaredGLMM(m)[1, "R2m"])
) |>
  dplyr::mutate(
    model      = factor(model, levels = names(decomp_models)),
    delta_R2m  = R2m - dplyr::lag(R2m, default = 0),
    component  = c("Intended change", "Observed compliance",
                   "Conscious detection", "Compliance × Detection")
  )

f2 <- ggplot(decomp_r2, aes(x = component, y = delta_R2m,
                              fill = component)) +
  geom_col(colour = "white", width = 0.65) +
  geom_text(aes(label = sprintf("+%.3f", delta_R2m)),
            vjust = -0.4, size = 3.5) +
  scale_fill_bcat(
                        guide = "none") +
  scale_x_discrete(limits = decomp_r2$component) +
  labs(
    title    = "Decomposition: unique R2m contribution per predictor",
    subtitle = "Incremental dR2m from each additional component",
    x = NULL, y = "dR2m (marginal R2)"
  ) +
  theme_bcat(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 9))

ggsave(file.path(belt_fig_dir, "S5_decomposition_R2.pdf"),
       f2, width = 8, height = 5, device = "pdf")
message("Saved: S5_decomposition_R2.pdf")


# ── F3: Within-person physio → arousal: hits vs misses ────────
if (nrow(within_comparison) == 0) {
  message("  [WARNING] within_comparison is empty — check term name matching. Skipping F3.")
} else {
  f3 <- within_comparison |>
    dplyr::mutate(
      condition = factor(condition, levels = c("Hits", "Misses"))
    ) |>
    ggplot(aes(x = condition, y = estimate,
               ymin = conf.low, ymax = conf.high,
               colour = condition)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.6) +
    geom_errorbar(width = 0.15, linewidth = 1.0) +
    geom_point(size = 4) +
    scale_colour_bcat(
                            guide = "none") +
    labs(
      title    = "Non-conscious compliance: physio -> arousal within-person",
      subtitle = "beta (dur_within -> Arousal) +/- 95% CI\nSimilar slopes = misattribution without awareness",
      x = "Detection outcome",
      y = "Within-person beta"
    ) +
    theme_bcat(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(belt_fig_dir, "S5_nonconscious_within.pdf"),
         f3, width = 5.5, height = 5, device = "pdf")
  message("Saved: S5_nonconscious_within.pdf")
}


# ── F4: Direction compliance rate by detection outcome ─────────
f4 <- dir_summary |>
  dplyr::mutate(
    Accuracy  = factor(Accuracy, levels = c(0, 1),
                        labels = c("Missed", "Detected")),
    se        = sqrt(pct_dir_correct * (100 - pct_dir_correct) /
                       n) / 100 * 100
  ) |>
  ggplot(aes(x = Accuracy, y = pct_dir_correct, fill = Accuracy)) +
  geom_hline(yintercept = 50, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  geom_col(width = 0.5, colour = "white") +
  geom_errorbar(aes(ymin = pct_dir_correct - se,
                    ymax = pct_dir_correct + se),
                width = 0.12, linewidth = 0.8) +
  annotate("text", x = 1, y = 52, hjust = 0.5, size = 3.5,
           label = sprintf("Binomial p = %.4f", binom_miss$p.value)) +
  scale_fill_bcat(
                        guide = "none") +
  scale_y_continuous(limits = c(0, 100),
                     labels = scales::percent_format(scale = 1)) +
  labs(
    title    = "Non-conscious compliance: direction correct rate",
    subtitle = "Dashed line = chance (50%)",
    x = "Detection outcome",
    y = "% trials breathing in correct direction"
  ) +
  theme_bcat(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(belt_fig_dir, "S5_direction_compliance.pdf"),
       f4, width = 5, height = 5, device = "pdf")
message("Saved: S5_direction_compliance.pdf")


# ── F5: Arousal × Change scatter split by session ─────────────
f5 <- regime_arousal$R1 |>
  dplyr::filter(!is.na(Arousal), !is.na(Change)) |>
  dplyr::mutate(
    Session = dplyr::if_else(run == 1, "Session 1", "Session 2"),
    Salience = dplyr::case_when(
      stringr::str_starts(tolower(condition), "high") ~ "High",
      stringr::str_starts(tolower(condition), "low")  ~ "Low",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(Salience)) |>
  ggplot(aes(x = Change, y = Arousal)) +
  geom_point(alpha = 0.06, size = 0.6, colour = "#2c7be5") +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2),
              colour = "#e63946", linewidth = 1.0, se = TRUE) +
  facet_grid(Session ~ Salience) +
  scale_x_continuous(
    breaks = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1\nFaster", "-0.5", "0\nNone", "0.5", "1\nSlower")
  ) +
  labs(
    title    = "H4A: Arousal by breathing change (both sessions, Regime 1)",
    subtitle = "Quadratic trend; each point = one trial",
    x = "Change (negative = faster)",
    y = "Arousal rating (1–6)"
  ) +
  theme_bcat(base_size = 12) +
  theme(plot.title   = element_text(face = "bold"),
        strip.text   = element_text(face = "bold"),
        panel.spacing = unit(1, "lines"))

ggsave(file.path(belt_fig_dir, "S5_H4A_arousal_change.pdf"),
       f5, width = 9, height = 7, device = "pdf")
message("Saved: S5_H4A_arousal_change.pdf")


# ============================================================
# 6. SALIENCE INDEPENDENCE CHECK
#
# Validates the theoretical claim that the High/Low Salience
# manipulation (abrupt vs. gradual onset) varies the detectability
# of the breathing change independently of its actual magnitude.
#
# If the claim holds, |dur_b1_vs_b4| should NOT differ by Salience
# after controlling for the intended change (delta). A small,
# non-significant Salience coefficient with Bayesian BF₀₁ > 3
# provides the empirical support for the independence claim.
#
# This analysis uses Regime 2 (trial-level QC applied) and both
# sessions. Results are saved to Results/ and referenced in the
# paper introduction.
# ============================================================
cat("\n\n========================================\n")
cat("6. SALIENCE INDEPENDENCE CHECK\n")
cat("Does Salience affect actual physiology?\n")
cat("Expected: NO — onset profile varies, not magnitude\n")
cat("========================================\n")

salience_check <- regime_arousal$R2 |>
  dplyr::filter(
    !is.na(dur_b1_vs_b4),
    !is.na(delta),
    delta != 0,       # change trials only; NoChange excluded
    !is.na(Salience)
  ) |>
  dplyr::mutate(
    abs_physio = abs(dur_b1_vs_b4),
    Salience   = factor(Salience, levels = c("Low", "High"))
  )

cat(sprintf("\nN trials for salience check: %d (%d participants)\n",
            nrow(salience_check), dplyr::n_distinct(salience_check$id)))

# Descriptive: mean |dur_b1_vs_b4| by Salience and Direction
sal_desc <- salience_check |>
  dplyr::group_by(Salience, Direction) |>
  dplyr::summarise(
    M_physio  = mean(abs_physio, na.rm = TRUE),
    SD_physio = sd(abs_physio,   na.rm = TRUE),
    n         = dplyr::n(),
    .groups   = "drop"
  )
cat("\nDescriptive: mean |physio change| by Salience × Direction:\n")
print(as.data.frame(sal_desc), digits = 3, row.names = FALSE)

# Primary model: does Salience predict physio after controlling for delta?
# delta is the intended change (same sign convention as dur_b1_vs_b4)
m_sal_main <- lmerTest::lmer(
  abs_physio ~ Salience + abs(delta) + (1 | id),
  data = salience_check, REML = FALSE
)

m_sal_null <- lmerTest::lmer(
  abs_physio ~ abs(delta) + (1 | id),
  data = salience_check, REML = FALSE
)

cat("\nMixed model: |dur_b1_vs_b4| ~ Salience + |delta| + (1|id)\n")
print(broom.mixed::tidy(m_sal_main, effects = "fixed", conf.int = TRUE),
      n = Inf)
cat(sprintf("  R²m = %.3f  R²c = %.3f\n",
            MuMIn::r.squaredGLMM(m_sal_main)[1],
            MuMIn::r.squaredGLMM(m_sal_main)[2]))

cat("\nLRT: does adding Salience improve fit over delta alone?\n")
lrt_sal <- anova(m_sal_null, m_sal_main)
print(lrt_sal)

# Bayesian null evidence for Salience effect
# Extract Salience t-value and df for BF approximation via BayesFactor
sal_cf  <- summary(m_sal_main)$coefficients
sal_row <- sal_cf["SalienceHigh", , drop = FALSE]
sal_t   <- sal_row[1, "t value"]
sal_df  <- sal_row[1, "df"]

# Point-null BF via t-distribution approximation (Rouder et al. 2009)
# For a simple check: r-scale = 0.707 (default JZS prior)
bf_sal <- tryCatch(
  BayesFactor::ttest.tstat(t = sal_t, n1 = round(sal_df) + 1,
                            nullInterval = NULL, complement = FALSE)$bf,
  error = function(e) NULL
)

if (!is.null(bf_sal)) {
  cat(sprintf("\nBayesian null evidence for Salience effect:\n"))
  cat(sprintf("  BF10 = %.3f  BF01 = %.3f\n", bf_sal, 1/bf_sal))
  cat(sprintf("  Interpretation: %s\n",
              dplyr::case_when(
                1/bf_sal >= 10 ~ "Strong evidence for null (independence confirmed)",
                1/bf_sal >= 3  ~ "Moderate evidence for null (independence supported)",
                1/bf_sal >= 1  ~ "Anecdotal evidence for null",
                TRUE           ~ "Evidence for Salience effect on physiology (check design)"
              )))
} else {
  cat("\n  [BF approximation unavailable — report LRT p-value only]\n")
}

# Robustness: also check Direction × Salience interaction
m_sal_int <- lmerTest::lmer(
  abs_physio ~ Salience * Direction + abs(delta) + (1 | id),
  data = salience_check |>
    dplyr::filter(Direction %in% c("Faster", "Slower")),
  REML = FALSE
)
cat("\nRobustness: Direction × Salience interaction on physiology:\n")
int_term <- grep("Salience.*Direction|Direction.*Salience",
                 rownames(summary(m_sal_int)$coefficients), value = TRUE)
if (length(int_term) > 0) {
  cat(sprintf("  %s\n", fmt_lmer(m_sal_int, int_term[1])))
} else {
  cat("  [No interaction term found]\n")
}

# Save results
sal_results <- tibble::tibble(
  analysis            = "Salience independence check",
  b_Salience          = sal_row[1, "Estimate"],
  se_Salience         = sal_row[1, "Std. Error"],
  t_Salience          = sal_t,
  df_Salience         = sal_df,
  p_Salience          = sal_row[1, "Pr(>|t|)"],
  BF01                = if (!is.null(bf_sal)) round(1/bf_sal, 3) else NA,
  lrt_p               = lrt_sal[["Pr(>Chisq)"]][2],
  n_trials            = nrow(salience_check),
  n_participants      = dplyr::n_distinct(salience_check$id),
  interpretation      = dplyr::case_when(
    sal_row[1, "Pr(>|t|)"] > .05 ~ "Salience does not predict physio magnitude (supports independence)",
    TRUE ~ "Salience predicts physio magnitude (check design)"
  )
)

readr::write_csv(sal_results,
                 file.path(RESULTS_DIR, "belt_salience_independence.csv"))
message("Saved: belt_salience_independence.csv")

# ============================================================
# 7. BELT SUMMARY FOR METHODS SECTION
# ============================================================
cat("\n\n========================================\n")
cat("BELT SUMMARY (methods reporting)\n")
cat("========================================\n")

belt_methods_summary <- qcSummary |>
  dplyr::filter(!is.na(run)) |>
  dplyr::group_by(belt_quality) |>
  dplyr::summarise(
    n_sessions = dplyr::n(),
    mean_pct_correct = mean(pct_correct_breaths, na.rm = TRUE),
    mean_r_delta     = mean(r_delta_change_all,  na.rm = TRUE),
    mean_pct_dir     = mean(pct_direction_correct, na.rm = TRUE),
    .groups = "drop"
  )

cat("Belt quality by session:\n")
print(as.data.frame(belt_methods_summary), digits = 3, row.names = FALSE)

regime_sizes <- tibble::tibble(
  Regime = names(regimes),
  n_participants = sapply(regime_arousal,
                          function(d) dplyr::n_distinct(d$id)),
  n_trials       = sapply(regime_arousal, nrow)
)

cat("\nRegime sample sizes:\n")
print(as.data.frame(regime_sizes), row.names = FALSE)

# Exclusions relative to R1
attr_r2b <- attr(regimes, "regime2b_exclude")
attr_r3  <- attr(regimes, "regime3_exclude")
cat(sprintf(
  "\nR2b excludes %d participants (r_delta<0.30 AND pct_correct<50)\n",
  length(attr_r2b)
))
cat(sprintf(
  "R3 additionally excludes %d participants (unusable belt)\n",
  length(setdiff(attr_r3, attr_r2b))
))

message("\nanalysis_belt.R complete.")
message("Figures saved to: ", belt_fig_dir)
