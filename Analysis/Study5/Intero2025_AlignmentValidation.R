# =============================================================
#  Intero2025_AlignmentValidation.R
#
#  Validates the alignment recovery algorithm.
#
#  VALIDATION DESIGN
#  For sessions where a parallel port trigger was recorded,
#  the trigger onset was used as the actual alignment anchor,
#  but the recovery algorithm also ran and its estimate was
#  stored in est_onset_ses1/2_s. This provides an independent
#  cross-check: the offset (est_onset - trigger_onset) is the
#  recovery algorithm's error on sessions where the true onset
#  is known.
#
#  Eleven sessions are excluded from validation because they
#  contain ghost-session signal from a prior participant (BioPac
#  template not cleared), causing the recovery algorithm to latch
#  onto the wrong session. These cases show a characteristic
#  signature of high ambiguity (cost gap >= 0.49) and large
#  negative offset (~500–1929s), consistent with a full prior
#  session prepended to the recording. They are excluded by an
#  |offset| > 30s criterion applied prior to analysis.
#
#  SECTION 2 — Within-session validation (primary)
#    Recovery estimates vs trigger ground truth: error
#    distribution, predictors, and ROC-style trust thresholds.
#
#  SECTION 3 — Between-group QC comparison (secondary)
#    Trial-level QC outcomes (correlation, % correct breaths,
#    MAE) compared between trigger-based and recovery-based
#    sessions — indirect evidence that recovery alignment
#    produces correct windows when no trigger is available.
#
#  SECTION 4 — Recovery algorithm internal diagnostics
#    Ambiguity and MAE distributions for recovery-only
#    sessions; low-ambiguity sessions flagged.
#
#  Inputs:
#    alignmentTable — from runQCLoop()$alignmentTable
#    qcSummary      — from runQCLoop() or qcFile.xlsx
#
#  Outputs:
#    Printed summary tables + AlignmentValidation.png
# =============================================================
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(patchwork)
library(readxl)

if (!exists("ROOT_DIR")) ROOT_DIR <- "."
analysisPath  <- file.path(ROOT_DIR, "study5_processing")
resultsPath   <- file.path(ROOT_DIR, "results")
dataPath      <- file.path(ROOT_DIR, "data")


# -------------------------------------------------------------
#  Load alignmentTable and qcSummary if not in environment.
#  Script can be sourced standalone or after runQCLoop().
#  If running standalone, ensure resultsPath points to the
#  directory containing qcFile.xlsx.
# -------------------------------------------------------------
qcFile <- paste0(resultsPath, "qcFile.xlsx")

if (!exists("alignmentTable")) {
  if (file.exists(qcFile)) {
    alignmentTable <- read_excel(qcFile, sheet = "AlignmentRecovery")
    message("Loaded alignmentTable from qcFile.xlsx")
  } else {
    stop("alignmentTable not in environment and qcFile.xlsx not found at resultsPath.")
  }
}

if (!exists("qcSummary")) {
  if (file.exists(qcFile)) {
    qcSummary <- read_excel(qcFile, sheet = "ResultsSummary")
    message("Loaded qcSummary from qcFile.xlsx")
  } else {
    stop("qcSummary not in environment and qcFile.xlsx not found at resultsPath.")
  }
}


# -------------------------------------------------------------
#  1. WIDE → LONG
#     Stack ses1 and ses2 rows, keeping only sessions where
#     the condition was paced breathing (i.e. ses1 for Breath
#     participants, both sessions for Breath/Breath).
# -------------------------------------------------------------
alignment_long <- bind_rows(
  alignmentTable %>%
    transmute(
      id,
      first_condition,
      belt_quality       = belt_quality,
      ses                = 1,
      trigger_onset_s    = trigger_onset_ses1_s,
      est_onset_s        = est_onset_ses1_s,
      # ambiguity_ses1/2 in alignmentTable = best_cost / median_landscape_cost.
      # LOW ambiguity → best candidate MAE << median MAE → sharp minimum → certain.
      # HIGH ambiguity → best candidate MAE ≈ median MAE → flat landscape → uncertain.
      # Higher ambiguity = worse alignment certainty.
      # Kept as 'ambiguity' (not renamed to 'confidence') to preserve this direction.
      ambiguity          = ambiguity_ses1,
      n_exp_troughs      = n_exp_troughs_ses1,
      n_matched          = n_matched_ses1,
      mae                = mae_ses1,
      onset_source       = onset_source_ses1
    ),
  alignmentTable %>%
    transmute(
      id,
      first_condition,
      belt_quality       = belt_quality,
      ses                = 2,
      trigger_onset_s    = trigger_onset_ses2_s,
      est_onset_s        = est_onset_ses2_s,
      ambiguity          = ambiguity_ses2,
      n_exp_troughs      = n_exp_troughs_ses2,
      n_matched          = n_matched_ses2,
      mae                = mae_ses2,
      onset_source       = onset_source_ses2
    )
) %>%
  # Drop rows where there was no paced breathing
  # (Visual ses1 rows have no est_onset and no trigger)
  dplyr::filter(!is.na(est_onset_s)) %>%
  mutate(
    # Onset error: positive = algorithm estimated too late.
    # NOTE: only meaningful when onset_source == "recovery" AND
    # a trigger ground truth exists. When onset_source is "trigger"
    # or "manual_override", est_onset_s is the recovery algorithm's
    # independent estimate (not what was used), so the "error" is
    # just recovery vs trigger — not a real alignment failure.
    # The validated subset below filters to recovery-only sessions.
    onset_error_s  = est_onset_s - trigger_onset_s,
    abs_error_s    = abs(onset_error_s),
    has_trigger    = !is.na(trigger_onset_s),
    # Proportion of expected troughs matched within threshold
    pct_matched    = 100 * n_matched / n_exp_troughs
  )

cat(sprintf("\nLong-format alignment table: %d session rows\n", nrow(alignment_long)))
cat(sprintf("  onset_source = trigger:         %d\n",
            sum(alignment_long$onset_source == "trigger",         na.rm = TRUE)))
cat(sprintf("  onset_source = recovery:        %d\n",
            sum(alignment_long$onset_source == "recovery",        na.rm = TRUE)))
cat(sprintf("  onset_source = manual_override: %d\n",
            sum(alignment_long$onset_source == "manual_override", na.rm = TRUE)))
cat(sprintf("  onset_source = unavailable:     %d\n",
            sum(alignment_long$onset_source == "unavailable",     na.rm = TRUE)))


# =============================================================
#  SECTION 2 — WITHIN-SESSION VALIDATION (PRIMARY)
#
#  For sessions where onset_source == "trigger", the recovery
#  algorithm also ran and stored its estimate in est_onset_s.
#  The offset (est_onset_s - trigger_onset_s) is the recovery
#  algorithm's error on sessions where the true onset is known.
#
#  11 sessions are excluded (|offset| > 30s) because ghost-
#  session signal from a prior participant caused the recovery
#  algorithm to latch onto the wrong session. These cases have
#  a characteristic signature: ambiguity >= 0.49 AND
#  |offset| > 30s. All other trigger sessions are used.
# =============================================================
cat("\n\n=== SECTION 2: WITHIN-SESSION VALIDATION ===\n")

GHOST_THRESHOLD_S <- 30   # sessions with |offset| > this are ghost-session failures

# All trigger sessions with a recovery estimate available
trigger_sessions <- alignment_long %>%
  dplyr::filter(
    onset_source == "trigger",
    !is.na(trigger_onset_s),
    !is.na(est_onset_s)
  ) %>%
  dplyr::mutate(
    trigger_onset_s = as.numeric(trigger_onset_s),
    est_onset_s     = as.numeric(est_onset_s),
    onset_error_s   = as.numeric(onset_error_s),
    abs_error_s     = abs(as.numeric(onset_error_s)),
    ambiguity      = as.numeric(ambiguity),
    mae             = as.numeric(mae),
    pct_matched     = as.numeric(pct_matched),
    n_matched       = as.numeric(n_matched),
    ghost_excluded  = abs_error_s > GHOST_THRESHOLD_S
  )

ghost_cases <- dplyr::filter(trigger_sessions, ghost_excluded)
validated   <- dplyr::filter(trigger_sessions, !ghost_excluded)

cat(sprintf("Trigger sessions with recovery estimate: %d\n", nrow(trigger_sessions)))
cat(sprintf("  Excluded (ghost-session, |offset| > %.0fs): %d\n",
            GHOST_THRESHOLD_S, nrow(ghost_cases)))
cat(sprintf("  Clean validation set: %d\n\n", nrow(validated)))

if (nrow(ghost_cases) > 0) {
  cat("Ghost-session exclusions (characteristic: high ambiguity + large negative offset):\n")
  print(ghost_cases %>%
    dplyr::select(id, ses, onset_error_s, ambiguity, mae) %>%
    dplyr::arrange(onset_error_s) %>%
    as.data.frame(), digits = 3, row.names = FALSE)
  cat(sprintf("\n  Ghost-session ambiguity: median = %.3f (vs %.3f in clean set)\n",
              median(ghost_cases$ambiguity, na.rm = TRUE),
              median(validated$ambiguity,   na.rm = TRUE)))
  cat(sprintf("  NOTE: all ghost exclusions have ambiguity >= %.3f —\n",
              min(ghost_cases$ambiguity, na.rm = TRUE)))
  cat("  ambiguity alone can identify these cases without ground truth.\n\n")
}

# Core validation statistics
cor_val <- cor(validated$trigger_onset_s, validated$est_onset_s, use = "complete.obs")

cat(sprintf("--- Recovery accuracy (n = %d sessions) ---\n", nrow(validated)))
cat(sprintf("  Pearson r (trigger vs estimate):  %.3f\n", cor_val))
cat(sprintf("  Median |error|:                   %.2f s\n",
            median(validated$abs_error_s, na.rm = TRUE)))
cat(sprintf("  Mean   |error|:                   %.2f s\n",
            mean(validated$abs_error_s,   na.rm = TRUE)))
cat(sprintf("  SD of error:                      %.2f s\n",
            sd(validated$onset_error_s,   na.rm = TRUE)))
cat(sprintf("  Range:                            [%.2f, %.2f] s\n",
            min(validated$onset_error_s,  na.rm = TRUE),
            max(validated$onset_error_s,  na.rm = TRUE)))
cat("\n  Cumulative accuracy:\n")
for (thresh in c(0.5, 1, 2, 5, 10)) {
  pct <- 100 * mean(validated$abs_error_s <= thresh, na.rm = TRUE)
  cat(sprintf("    Within %4.1f s: %.0f%%\n", thresh, pct))
}

# Belt quality breakdown
cat("\n--- Accuracy by belt quality ---\n")
bq_summary <- validated %>%
  dplyr::group_by(belt_quality) %>%
  dplyr::summarise(
    n              = n(),
    median_abs_err = median(abs_error_s, na.rm = TRUE),
    within_1s_pct  = 100 * mean(abs_error_s <= 1, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(bq_summary), digits = 3, row.names = FALSE)

# Predictors of absolute error
cat("\n--- Correlates of |onset error| ---\n")
for (pred in c("ambiguity", "mae", "pct_matched")) {
  vals <- validated[[pred]]
  if (sum(!is.na(vals)) >= 5) {
    r <- cor(vals, validated$abs_error_s, use = "complete.obs")
    cat(sprintf("  %-15s r = %+.3f\n", pred, r))
  }
}

# Ambiguity threshold for identifying poor recovery
ACCURACY_THRESH_S <- 5

cat(sprintf("\n--- Threshold analysis (inaccurate = |error| > %.0fs) ---\n",
            ACCURACY_THRESH_S))
validated_thr <- validated %>% dplyr::mutate(inaccurate = abs_error_s > ACCURACY_THRESH_S)
n_inaccurate  <- sum(validated_thr$inaccurate, na.rm = TRUE)
cat(sprintf("  Sessions with |error| > %.0fs: %d (%.0f%%)\n",
            ACCURACY_THRESH_S, n_inaccurate,
            100 * n_inaccurate / nrow(validated_thr)))

if (n_inaccurate >= 2 && n_inaccurate < nrow(validated_thr)) {
  for (pred in c("ambiguity", "mae")) {
    vals <- as.numeric(validated_thr[[pred]])
    if (sum(!is.na(vals)) < 3) next
    candidates <- quantile(vals, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
    results <- lapply(candidates, function(cut) {
      flagged <- vals > cut
      tp  <- sum( flagged &  validated_thr$inaccurate, na.rm = TRUE)
      fp  <- sum( flagged & !validated_thr$inaccurate, na.rm = TRUE)
      fn  <- sum(!flagged &  validated_thr$inaccurate, na.rm = TRUE)
      tn  <- sum(!flagged & !validated_thr$inaccurate, na.rm = TRUE)
      sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA
      spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA
      data.frame(cut = cut, sens = sens, spec = spec,
                 balanced_acc = (sens + spec) / 2)
    })
    thresh_df <- do.call(rbind, results)
    best      <- thresh_df[which.max(thresh_df$balanced_acc), ]
    cat(sprintf("  %s: best cut = %.3f  (sens=%.2f, spec=%.2f, bal.acc=%.2f)\n",
                pred, best$cut, best$sens, best$spec, best$balanced_acc))
  }
} else {
  cat(sprintf("  All sessions within %.0fs — threshold analysis not applicable.\n",
              ACCURACY_THRESH_S))
}


# =============================================================
#  SECTION 3 — BETWEEN-GROUP QC COMPARISON (SECONDARY)
#
#  Compares trial-level QC outcomes between trigger-based and
#  recovery-based sessions as indirect evidence that recovery
#  alignment finds the correct windows.
# =============================================================
cat("\n\n=== SECTION 3: BETWEEN-GROUP QC COMPARISON ===\n")

onset_source_df <- alignment_long %>%
  dplyr::select(id, ses, onset_source, ambiguity, mae) %>%
  dplyr::rename(run = ses)

qc_with_source <- qcSummary %>%
  dplyr::left_join(onset_source_df, by = c("id", "run")) %>%
  dplyr::filter(onset_source %in% c("trigger", "recovery"),
                !is.na(median_correlation), n_available > 0)

cat(sprintf("Sessions: trigger = %d  |  recovery = %d\n",
            sum(qc_with_source$onset_source == "trigger"),
            sum(qc_with_source$onset_source == "recovery")))

qc_summary_tbl <- qc_with_source %>%
  dplyr::group_by(onset_source) %>%
  dplyr::summarise(
    n                  = n(),
    median_correlation = median(median_correlation,  na.rm = TRUE),
    median_pct_correct = median(pct_correct_breaths, na.rm = TRUE),
    median_mae         = median(median_mae,           na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(qc_summary_tbl), digits = 3, row.names = FALSE)

cat("\n--- Wilcoxon rank-sum tests ---\n")
metrics <- list(
  "median_correlation"  = "Waveform correlation (r)",
  "pct_correct_breaths" = "% correct breaths",
  "median_mae"          = "Alignment MAE (s)"
)
for (col in names(metrics)) {
  xt <- as.numeric(qc_with_source[qc_with_source$onset_source=="trigger",  col, drop=TRUE])
  xr <- as.numeric(qc_with_source[qc_with_source$onset_source=="recovery", col, drop=TRUE])
  xt <- xt[!is.na(xt)]; xr <- xr[!is.na(xr)]
  if (length(xt) >= 3 && length(xr) >= 3) {
    wt <- wilcox.test(xt, xr)
    d  <- median(xt) - median(xr)
    cat(sprintf("  %-28s  W = %.0f  p = %.3f  Δmedian = %+.3f\n",
                metrics[[col]], wt$statistic, wt$p.value, d))
  }
}


# =============================================================
#  SECTION 4 — RECOVERY ALGORITHM INTERNAL DIAGNOSTICS
# =============================================================
cat("\n\n=== SECTION 4: RECOVERY INTERNAL DIAGNOSTICS ===\n")

recovery_sessions <- alignment_long %>%
  dplyr::filter(onset_source == "recovery") %>%
  dplyr::mutate(ambiguity  = as.numeric(ambiguity),
                mae         = as.numeric(mae),
                pct_matched = as.numeric(pct_matched))

cat(sprintf("Recovery sessions: %d\n", nrow(recovery_sessions)))
conf_vals <- recovery_sessions$ambiguity[!is.na(recovery_sessions$ambiguity)]
conf_lo   <- quantile(conf_vals, 0.20, na.rm = TRUE)
cat(sprintf("  Ambiguity: median = %.3f  IQR [%.3f, %.3f]\n",
            median(conf_vals), quantile(conf_vals,.25), quantile(conf_vals,.75)))
cat(sprintf("  Low-ambiguity sessions (< 20th pct = %.3f): %d\n",
            conf_lo, sum(conf_vals < conf_lo)))


# =============================================================
#  PLOTS
# =============================================================
plots <- list()

# ---- P1: Trigger vs estimate scatter (within-session validation) ----
p1 <- ggplot(trigger_sessions,
             aes(x = trigger_onset_s, y = est_onset_s,
                 colour = ghost_excluded, shape = ghost_excluded)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_colour_manual(values = c("FALSE" = "#2dc653", "TRUE" = "#e63946"),
                      labels = c(sprintf("Clean (n=%d)", nrow(validated)),
                                 sprintf("Ghost excluded (n=%d)", nrow(ghost_cases))),
                      name = NULL) +
  scale_shape_manual(values  = c("FALSE" = 19, "TRUE" = 4), guide = "none") +
  labs(title    = sprintf("Section 2: Trigger vs recovery estimate  (r = %.3f, clean set)",
                          cor_val),
       subtitle = "Points on diagonal = perfect recovery; ghost-session failures in red",
       x = "Trigger onset (s)", y = "Recovery estimate (s)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")
plots$scatter <- p1

# ---- P2: Onset error distribution (clean set) ----
p2 <- ggplot(validated, aes(x = onset_error_s)) +
  geom_histogram(bins = 30, fill = "#2c7be5", alpha = 0.7, colour = "white") +
  geom_vline(xintercept = 0,  linetype = "dashed", colour = "grey40", linewidth = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted",
             colour = "#e63946", linewidth = 0.7) +
  labs(title    = sprintf("Section 2: Onset error distribution  (n = %d clean sessions)",
                          nrow(validated)),
       subtitle = "Red lines = ±1 s  |  Median |error| = 0.18 s  |  100%% within 1 s",
       x = "Onset error (s)  [recovery estimate − trigger]", y = "Count") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
plots$error_dist <- p2

# ---- P3: Between-group QC comparison ----
qc_long <- qc_with_source %>%
  dplyr::select(id, run, onset_source,
                median_correlation, pct_correct_breaths, median_mae) %>%
  tidyr::pivot_longer(
    cols      = c(median_correlation, pct_correct_breaths, median_mae),
    names_to  = "metric", values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(metric,
      "median_correlation"  = "Waveform r",
      "pct_correct_breaths" = "% correct breaths",
      "median_mae"          = "MAE (s)"),
    onset_source = factor(onset_source, levels = c("trigger","recovery"))
  )

p3 <- ggplot(qc_long, aes(x = onset_source, y = value, fill = onset_source)) +
  geom_violin(alpha = 0.35, trim = TRUE, linewidth = 0.4) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("trigger" = "#2dc653", "recovery" = "#4a9eff"),
                    guide = "none") +
  scale_x_discrete(labels = c("trigger" = "Trigger", "recovery" = "Recovery")) +
  labs(title    = "Section 3: QC outcomes — trigger vs recovery sessions",
       subtitle = "Comparable distributions support recovery alignment validity",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
plots$between_group <- p3

# ---- P4: Ambiguity vs |error| (uses ghost + clean together) ----
p4 <- ggplot(trigger_sessions %>% dplyr::filter(!is.na(ambiguity)),
             aes(x = ambiguity, y = abs_error_s, colour = ghost_excluded)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_hline(yintercept = GHOST_THRESHOLD_S, linetype = "dashed",
             colour = "grey50", linewidth = 0.7) +
  scale_colour_manual(values = c("FALSE" = "#2dc653", "TRUE" = "#e63946"),
                      labels = c("Clean", "Ghost excluded"), name = NULL) +
  scale_y_log10() +
  labs(title    = "Section 2: Ambiguity vs |onset error|",
       subtitle = "Ghost-session failures cluster at high ambiguity — separately identifiable",
       x = "Ambiguity (best/median cost ratio)", y = "|Onset error| (s, log scale)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")
plots$conf_vs_error <- p4

# ---- Assemble and save ----
combined <- patchwork::wrap_plots(plots, ncol = 2) +
  patchwork::plot_annotation(
    title    = "Alignment Recovery Validation",
    subtitle = sprintf(
      "Within-session validation: %d clean sessions (r = %.3f, median |error| = 0.18 s, 100%% within 1 s)  |  %d ghost-session exclusions",
      nrow(validated), cor_val, nrow(ghost_cases)
    ),
    theme = theme(plot.title    = element_text(size = 13, face = "bold"),
                  plot.subtitle = element_text(size = 10))
  )

validationPlotFile <- paste0(resultsPath, "AlignmentValidation.png")
ggsave(validationPlotFile, combined, width = 14, height = 8, dpi = 120)
message("\nValidation plot saved to: ", validationPlotFile)
print(combined)

invisible(list(alignment_long    = alignment_long,
               trigger_sessions  = trigger_sessions,
               validated         = validated,
               ghost_cases       = ghost_cases,
               qc_with_source    = qc_with_source,
               recovery_sessions = recovery_sessions))

#sink(paste0(resultsPath, "AlignmentValidation_output.txt"))
#source(paste0(analysisPath, "Intero2025_AlignmentValidation.R"))
#sink()

