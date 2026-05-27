# ============================================================
# MainAnalysis.R
# BCAT Five-Study Paper — Master Script
#
# Sets all paths once. Every analysis file reads BASE_DIR,
# DATA_DIR, ANALYSIS_DIR, RESULTS_DIR, FIG_DIR, MODEL_DIR
# from this environment — no paths are set anywhere else.
#
# Run order:
#   1. DataCleaning/study{N}_clean.R       (data preparation, run via run_all_cleaning.R)
#   2. Study5/Intero2025_PrepQualtrics.R   (Study 5 questionnaires)
#   3. This script                         (sources all analysis files in order)
#
# UTILITIES (sourced first, no output):
#   utils.R                               — helper functions and load_all_data()
#   theme_bcat.R                          — shared ggplot theme and save_bcat_fig()
#   meta_analysis.R                       — meta-analytic functions (run_h4b_meta() etc.)
#
# PRIMARY RESULTS (main text order):
#   analysis_val_detection.R              — Task Validation: Change² detection models
#                                           Creates: table_detection_change2.csv
#   analysis_val_pilot_studies.R          — Supplement S1.1–S1.4: pilot comparison,
#                                           staircase convergence, salience orthogonality
#                                           Creates: study1_task_comparison.csv,
#                                                    table_study3_salience_accuracy.csv,
#                                                    table_staircase_convergence.csv,
#                                                    study2_convergence_elbow.csv,
#                                                    study2_convergence_elbow_summary.csv
#   analysis_val_thresholds.R             — Supplement S1.5: thresholds, d', ICC (H1, H2, H5, H6)
#                                           Creates: table_validation.csv,
#                                                    table_reliability.csv
#   analysis_arousal.R                    — Awareness Gates Arousal Transfer (H4A, H4B, H4C)
#                                           Creates: table_arousal.csv,
#                                                    table_h4c_exploratory.csv,
#                                                    meta_h4b_sensitivity_noS3.csv
#   analysis_miss_baseline.R              — Supplement S2.1: Bayesian null, misses vs. baseline
#                                           Creates: miss_baseline_bf.csv
#   test_block_arousal.R                  — Supplement S2.2: test block mediation
#                                           Creates: test_block_all_models.csv,
#                                                    frequentist_mediation_by_group.csv,
#                                                    brms_indirect_by_group.csv,
#                                                    brms_arousal_indirect_draws.csv,
#                                                    freq_vs_brms_key_terms.csv,
#                                                    test_block_study4_stratified.csv,
#                                                    test_block_study5_sessions.csv
#   analysis_study3_attraction_mediation.R — Study 3 misattribution chain
#                                           Creates: study3_mediation_primary.csv,
#                                                    study3_mediation_m1a_subset.csv,
#                                                    study3_mediation_m1b_analytic.csv,
#                                                    study3_mediation_bayesian_m0.csv,
#                                                    study3_mediation_bayesian_m1.csv,
#                                                    study3_maia_moderation_path_a.csv,
#                                                    study3_maia_protection_acme.csv,
#                                                    study3_maia_bayesian_moderation.csv
#   analysis_maia.R                       — MAIA: sensibility vs. sensitivity (H3A, H3B)
#                                           Creates: table_maia.csv,
#                                                    table_maia_gating_moderation.csv,
#                                                    table_threshold_descriptives.csv,
#                                                    study4_bh_correction.csv,
#                                                    study5_bh_correction.csv
#
# SUPPLEMENTARY (supplement order):
#   analysis_tce.R                        — S3: TCE sensitivity analyses (regime, matched, prior);
#                                           Study 5 completer check
#                                           Creates: tce_primary_results.csv,
#                                                    barrett_tce_bayes.csv,
#                                                    tce_sensitivity_regime.csv,
#                                                    tce_sensitivity_matched.csv,
#                                                    tce_sensitivity_prior.csv,
#                                                    tce_sensitivity_consolidated.csv,
#                                                    s5_completer_check.csv
#   analysis_belt.R                       — S4: belt physio and compliance
#                                           Creates: table_belt_compliance.csv,
#                                                    table_belt_physio_arousal.csv,
#                                                    belt_direction_compliance_misses.csv,
#                                                    belt_regime_key_terms.csv,
#                                                    belt_salience_independence.csv
#   belt_salience_followup.R              — S4: belt salience follow-up (exploratory)
#                                           Creates: table_belt_salience_followup.csv
#   analysis_s4_entrainment.R             — S4: belt-pacer entrainment -> detection advantage
#                                           Creates: s4_entrainment_results.csv
#   analysis_individual_differences.R     — S5: individual differences (Studies 4 and 5)
#                                           Creates: id_features_{study}.csv,
#                                                    id_correlations_study4.csv,
#                                                    id_correlations_study5.csv,
#                                                    id_replication_study4_5.csv
#   analysis_study5_exploratory.R         — S5: exploratory Study 5 MAIA and threshold breakdown
#                                           Creates: s5e_h3b_threshold_breakdown.csv,
#                                                    s5e_maia_subscale_cors.csv,
#                                                    s5e_alexithymia_dissociation.csv,
#                                                    s5e_group_awareness.csv,
#                                                    s5e_training_contrast.csv,
#                                                    s5e_awareness_change.csv
#   analysis_s7_maia_selfesteem.R         — S7: MAIA-confidence controlling for self-esteem
#                                           Creates: s7_partial_correlations.csv,
#                                                    s7_threshold_partial_correlations.csv,
#                                                    s7_multilevel_confidence_results.csv,
#                                                    s7_multilevel_threshold_results.csv
#   analysis_hbd.R                        — S6: heartbeat detection, cross-modal dissociation
#                                           Creates: hbd_missingness.csv,
#                                                    hbd3_crossmodal.csv,
#                                                    hbd5_sensitivity.csv,
#                                                    hbd_summary.csv
#
# FIGURES (sourced after all analyses):
#   fig_staircase.R                       — Figure S1 (staircase convergence)
#   fig_accuracy.R                        — Figure 1 (psychometric functions)
#   fig_arousal.R                         — Figure 2 (arousal gating)
#   fig_regime_comparison.R               — Figure S3 (TCE regime sensitivity)
#
# TABLES (sourced at end):
#   Build_Main_Tables.R                   — assembles BCAT_Main_Tables.docx from CSVs
#   Build_Reliability_Tables.R            — assembles Scale_Reliability_Tables.docx
#   Build_Supplementary_Tables.R          — assembles BCAT_Supplementary_Tables.docx
#
# STANDALONE (not sourced here — run separately):
#   ScaleReliability/Study{N}_PrepScales.R — scale scoring per study
#   DataCleaning/run_all_cleaning.R       — runs all study{N}_clean.R scripts
#   Study5/                               — full physio pipeline (see Study5/study5_processing_README.md)
#
# Shared objects available to all sourced files:
#   Data:    s1l, s1s, s2l, s2s, s3l, s3s, s4l, s4s, s4t,
#            s5l, s5s, s5t, s5_long_breath, s1l_A_sorted
#   Paths:   BASE_DIR, DATA_DIR, ANALYSIS_DIR, RESULTS_DIR,
#            FIG_DIR, TABLE_DIR, MODEL_DIR, RDS_DIR
#   Helpers: partial_r_from_t(), add_odds_ratios(), add_partial_r(),
#            extract_change2_results(), compute_test_dprime_3afc(),
#            make_threshold_long(), compute_retest_icc(), %||%
#            apply_bh_correction() — all defined in utils.R
# ============================================================

# ── Packages ──────────────────────────────────────────────────
packages <- c(
  # Data I/O
  "tidyverse", "readxl",
  # Mixed models
  "lme4", "lmerTest",
  # Model utilities
  "emmeans", "broom.mixed", "MuMIn", "irr",
  # Statistical tests
  "BayesFactor", "ppcor", "car",
  # Bayesian
  "brms", "posterior",
  # Meta-analysis
  "metafor",
  # Mediation
  "mediation",
  # Visualization
  "patchwork", "ggeffects", "scales",
  # Table building (Build_*.R files sourced below)
  "flextable", "officer"
)
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
options(readr.show_col_types = FALSE,
        dplyr.summarise.inform = FALSE)

for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}


# ============================================================
# Paths — set once here; all sourced files use these objects
# ============================================================
BASE_DIR     <- "C:/Users/norma/Desktop/test6/" # <-- insert your base path here!!!
DATA_DIR     <- file.path(BASE_DIR, "Data")
ANALYSIS_DIR <- file.path(BASE_DIR, "Analysis")
RESULTS_DIR  <- file.path(BASE_DIR, "Results")
FIG_DIR      <- file.path(BASE_DIR, "Figures")
TABLE_DIR      <- file.path(BASE_DIR, "Tables")
MODEL_DIR    <- file.path(RESULTS_DIR, "Models")
RDS_DIR <- file.path(RESULTS_DIR, "Models", "TestRDS")

for (.d in c(RESULTS_DIR, FIG_DIR, TABLE_DIR, MODEL_DIR, RDS_DIR)) {
  dir.create(.d, showWarnings = FALSE, recursive = TRUE)
}
rm(.d)

# ── Utility scripts ───────────────────────────────────────────
# utils.R:        data loading, effect size helpers, reshape helpers
# theme_bcat.R:   shared ggplot theme
# meta_analysis.R: meta-analytic pooling helpers
source(file.path(ANALYSIS_DIR, "utils.R"))
source(file.path(ANALYSIS_DIR, "theme_bcat.R"))
source(file.path(ANALYSIS_DIR, "meta_analysis.R"))


# ── Data loading and standardisation ──────────────────────────
message("Loading and standardising data...")

d <- load_all_data()

for (study in 1:5) {
  key_l <- paste0("s", study, "_long")
  key_s <- paste0("s", study, "_summary")
  if (key_l %in% names(d))
    d[[key_l]] <- standardise_study_data(d[[key_l]], study)
  if (key_s %in% names(d))
    d[[key_s]] <- standardise_study_data(d[[key_s]], study) |>
      standardise_maia(study_label = paste("Study", study))
}

s1l <- d$s1_long;  s1s <- d$s1_summary
s2l <- d$s2_long;  s2s <- d$s2_summary
s3l <- d$s3_long;  s3s <- d$s3_summary
s4l <- d$s4_long;  s4s <- d$s4_summary;  s4t <- d$s4_test
s5l <- d$s5_long;  s5s <- d$s5_summary;  s5t <- d$s5_test

# Study 5 belt QC (used by analysis_belt.R, belt_salience_followup.R,
#                   analysis_s6_entrainment.R)
qcFull    <- d$s5_qcFull
qcSummary <- d$s5_qcSummary

# Study 5 HBD (used by analysis_hbd.R)
s5_hbd           <- d$s5_hbd
s5_hbd_intervals <- d$s5_hbd_intervals

# Parse Salience and Direction from Condition strings in test files
s4t <- s4t |>
  dplyr::mutate(
    Salience  = dplyr::if_else(stringr::str_starts(Condition,     "high"), "High", "Low"),
    Direction = dplyr::if_else(stringr::str_ends(Condition,       "Acc"),  "Faster", "Slower")
  )
s5t <- s5t |>
  dplyr::mutate(
    Salience  = dplyr::if_else(stringr::str_starts(taskCondition, "high"), "High", "Low"),
    Direction = dplyr::if_else(stringr::str_ends(taskCondition,   "Acc"),  "Faster", "Slower")
  )

# Study 5: session-mean Confidence and Awareness across sessions
s5s <- s5s |>
  dplyr::mutate(
    mean_Confidence = rowMeans(
      dplyr::across(c(mean_Confidence_ses1, mean_Confidence_ses2)), na.rm = TRUE),
    Awareness = rowMeans(
      dplyr::across(c(Awareness_ses1, Awareness_ses2)), na.rm = TRUE)
  )

# Convenience subsets used across multiple analysis files
s5_long_breath <- dplyr::filter(s5l, Condition == "breath")

# Pre-sorted Study 1A for MRC threshold computation
s1l_A_sorted <- s1l |>
  dplyr::filter(Group == "TaskA") |>
  dplyr::arrange(id, TrialNum)

message("Data ready.")


# ── PRIMARY RESULTS ───────────────────────────────────────────

# Task Validation: Change² dose-response across all studies (main text Table 2)
# Creates: table_detection_change2.csv
# Creates env objects: s{N}_det — used by analysis_val_thresholds.R
source(file.path(ANALYSIS_DIR, "analysis_val_detection.R"))

# Task Validation: pilot study comparisons and staircase convergence (Supplement S1.1–S1.4)
# Creates: study1_task_comparison.csv, table_study3_salience_accuracy.csv,
#          table_staircase_convergence.csv, study2_convergence_elbow.csv,
#          study2_convergence_elbow_summary.csv
source(file.path(ANALYSIS_DIR, "analysis_val_pilot_studies.R"))

# Task Validation: thresholds, d', ICC (Supplement S1.4)
# NOTE: requires analysis_val_detection.R to have run first
#       (uses s{N}_det objects for table_validation.csv assembly)
# Creates: table_validation.csv, table_reliability.csv
# Creates env objects: s{N}_thresh, s{N}_conf, s{N}_maia — used by analysis_maia.R
source(file.path(ANALYSIS_DIR, "analysis_val_thresholds.R"))

# Test block detection accuracy: salience, direction, group, session
# Generates table_test_dprime by group/direction (replaces table_test_dprime.csv)
# Creates: test_block_accuracy_descriptives.csv, test_block_dprime_by_group.csv,
#          test_block_dprime_summary.csv, test_block_accuracy_glmm.csv,
#          test_block_dprime_lmm.csv, test_block_study4_accuracy.csv,
#          test_block_study5_accuracy.csv, test_block_direction_asymmetry.csv
source(file.path(ANALYSIS_DIR, "test_block_accuracy.R"))

# Awareness Gates Arousal Transfer (H4A, H4B, H4C)
# Creates: table_arousal.csv, meta_h4b_sensitivity_noS3.csv,
#          table_h4c_exploratory.csv
# Creates env objects: s{N}_arousal, meta_h4b — used by analysis_maia.R
source(file.path(ANALYSIS_DIR, "analysis_arousal.R"))

# Bayesian null test: missed-change trials vs. no-change baseline
# BF01 values reported in main text (Studies 1A, 2, 4, 5)
# Creates: miss_baseline_bf.csv
source(file.path(ANALYSIS_DIR, "analysis_miss_baseline.R"))

# Supplement S2.2: test block mediation (salience -> detection -> arousal/confidence)
# Creates: test_block_all_models.csv, frequentist_mediation_by_group.csv,
#          brms_indirect_by_group.csv, brms_arousal_indirect_draws.csv,
#          freq_vs_brms_key_terms.csv, test_block_study4_stratified.csv,
#          test_block_study5_sessions.csv
source(file.path(ANALYSIS_DIR, "test_block_arousal.R"))

# Study 3 misattribution chain: Change -> Arousal -> Attraction
# Creates: study3_mediation_primary.csv, study3_mediation_m1a_subset.csv,
#          study3_mediation_m1b_analytic.csv, study3_mediation_bayesian_m0.csv,
#          study3_mediation_bayesian_m1.csv, study3_maia_moderation_path_a.csv,
#          study3_maia_protection_acme.csv, study3_maia_bayesian_moderation.csv
source(file.path(ANALYSIS_DIR, "analysis_study3_attraction_mediation.R"))

# MAIA: sensibility vs. sensitivity (H3A, H3B)
# NOTE: requires analysis_val_thresholds.R and analysis_arousal.R to have run first
#       (uses s{N}_maia, s{N}_thresh, s{N}_conf, s{N}_arousal, meta_h4b objects)
# Creates: table_maia.csv, table_maia_gating_moderation.csv,
#          table_threshold_descriptives.csv, study4_bh_correction.csv,
#          study5_bh_correction.csv
source(file.path(ANALYSIS_DIR, "analysis_maia.R"))


# ── SUPPLEMENTARY ─────────────────────────────────────────────

# S3: TCE sensitivity analyses (regime, matched, prior); Study 5 completer check
# Creates: tce_primary_results.csv, barrett_tce_bayes.csv,
#          tce_sensitivity_regime.csv, tce_sensitivity_matched.csv,
#          tce_sensitivity_prior.csv, tce_sensitivity_consolidated.csv,
#          s5_completer_check.csv
source(file.path(ANALYSIS_DIR, "analysis_tce.R"))

# S4: Belt physio, compliance, salience independence, regime comparison
# Creates: table_belt_compliance.csv, table_belt_physio_arousal.csv,
#          belt_direction_compliance_misses.csv, belt_regime_key_terms.csv,
#          belt_salience_independence.csv
source(file.path(ANALYSIS_DIR, "analysis_belt.R"))
source(file.path(ANALYSIS_DIR, "belt_salience_followup.R"))
source(file.path(ANALYSIS_DIR, "analysis_s4_entrainment.R"))

# S5: Individual differences (Studies 4 and 5 replication)
source(file.path(ANALYSIS_DIR, "analysis_individual_differences.R"))

# Exploratory Study 5 MAIA and threshold breakdown
source(file.path(ANALYSIS_DIR, "analysis_study5_exploratory.R"))

# S7: Controls for Global Self Esteem on MAIA in Study 5
source(file.path(ANALYSIS_DIR, "analysis_s7_maia_selfesteem.R")) 

# HBD cardiac interoception (exploratory)
source(file.path(ANALYSIS_DIR, "analysis_hbd.R"))

# Figures
source(file.path(ANALYSIS_DIR, "fig_staircase.R"))
source(file.path(ANALYSIS_DIR, "fig_accuracy.R"))
source(file.path(ANALYSIS_DIR, "fig_arousal.R"))
source(file.path(ANALYSIS_DIR, "fig_regime_comparison.R"))

# Tables
source(file.path(ANALYSIS_DIR, "Build_Main_Tables.R"))
source(file.path(ANALYSIS_DIR, "Build_Reliability_Tables.R"))
source(file.path(ANALYSIS_DIR, "Build_Supplementary_Tables.R"))

# ── Zip all results ───────────────────────────────────────────
result_files <- list.files(RESULTS_DIR, pattern = "\\.csv$", full.names = TRUE)
zip(zipfile = file.path(RESULTS_DIR, "all_results.zip"), files = result_files)
message("Zipped ", length(result_files), " CSV files to all_results.zip")

message("\nMainAnalysis.R complete.")
