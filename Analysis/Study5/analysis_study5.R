# =============================================================
#  analysis_study5.R
#  Study 5 — In-person BCAT with respiratory belt
#
#  Preregistered hypotheses H1–H6 with sensitivity analysis
#  across four physio exclusion regimes.
#
#  INPUT FILES
#  ─────────────────────────────────────────────────────────
#  study5_long.csv     — trial-level staircase data
#  study5_summary.csv  — participant-level thresholds, d',
#                        MAIA, and questionnaire scores
#  study5_test.csv     — test-phase trial data (H5)
#  qcFile.xlsx         — belt QC output (FullResults +
#                        ResultsSummary)
#
#  REGIME STRUCTURE
#  ─────────────────────────────────────────────────────────
#  Regime 1 (primary): full physio sample, no exclusions
#  Regime 2: trial-level conjunction rule (≥2 of 3 flags)
#  Regime 2b: + session-level manipulation check exclusions
#  Regime 3: + unusable-belt participant exclusions
#
#  HYPOTHESIS STRUCTURE
#  ─────────────────────────────────────────────────────────
#  H1, H2, H3, H5, H6 — behavioural only; run once on full
#    sample (physio regime does not apply)
#  H4A–C — arousal × respiratory change; run across all
#    four regimes as sensitivity analysis
#
#  SIGN CONVENTION
#  ─────────────────────────────────────────────────────────
#  Change (from study5_long.csv): negative = Faster (higher
#    breathing frequency); positive = Slower.
#  The preregistration predicts faster breathing → higher
#    arousal, so H4A expects a negative β on Change.
#  Change is used directly as the predictor throughout
#  (no re-centring needed; 0 = NoChange trials).
#
#  MULTIPLE COMPARISONS
#  ─────────────────────────────────────────────────────────
#  Benjamini–Hochberg correction applied in two families:
#    Family 1 — main effects: H1, H2, H4A, H4C-interaction,
#               H5 (one critical p-value each)
#    Family 2 — interaction terms: H4B, H4C
#  H3 uses Bayesian BF for the null; excluded from BH.
#  H6 is a reliability estimate; excluded from BH.
# =============================================================

# =============================================================
#  SET UP
# =============================================================
packages <- c(
  "tidyverse", "readxl", "lme4", "lmerTest",
  "broom.mixed", "BayesFactor", "patchwork", "viridis"
)
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
options(readr.show_col_types = FALSE)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

# Paths — update to match your environment
if (!exists("ROOT_DIR")) ROOT_DIR <- "."
analysisPath  <- file.path(ROOT_DIR, "study5_processing")
resultsPath   <- file.path(ROOT_DIR, "results")
dataPath      <- file.path(ROOT_DIR, "data")



# =============================================================
#  1. LOAD RAW DATA
# =============================================================

# ---- 1a. Belt QC ----
qcFile    <- file.path(resultsPath, "qcFile.xlsx")
qcFull    <- readxl::read_excel(qcFile, sheet = "FullResults")
qcSummary <- readxl::read_excel(qcFile, sheet = "ResultsSummary")

# ---- 1b. Behavioural files ----
long    <- readr::read_csv(file.path(resultsPath, "study5_long.csv"))
summary_data <- readr::read_csv(file.path(dataPath, "study5_summary.csv"))
test    <- readr::read_csv(file.path(resultsPath, "study5_test.csv"))


# =============================================================
#  2. DATA CLEANING AND PREPARATION
# =============================================================

# ---- 2a. Fix Condition NAs in long ----
# Condition was not populated for ses2 in Visual-first participants.
# Reconstruct from Group: Group=="breath" → ses1=Breath, ses2=Visual
#                         Group=="visual" → ses1=Visual, ses2=Breath
long <- long %>%
  dplyr::left_join(
    summary_data %>% dplyr::select(id, Group),
    by = "id"
  ) %>%
  dplyr::mutate(
    Condition = dplyr::case_when(
      !is.na(Condition)                       ~ Condition,
      Group == "breath" & ses == "ses1"       ~ "breath",
      Group == "breath" & ses == "ses2"       ~ "breath",
      Group == "visual" & ses == "ses1"       ~ "visual",
      Group == "visual" & ses == "ses2"       ~ "breath",
      TRUE                                    ~ NA_character_
    ),
    # Recode ses → run integer to match qcFile
    run = dplyr::if_else(ses == "ses1", 1L, 2L)
  )

# ---- 2b. Fix test file: parse Salience and Direction from taskCondition ----
# Salience and Direction columns are entirely NA in test file.
# Reconstruct from taskCondition (e.g. "highSalienceAcc" → High / Faster)
test <- test %>%
  dplyr::mutate(
    Salience  = dplyr::case_when(
      stringr::str_starts(taskCondition, "high") ~ "High",
      stringr::str_starts(taskCondition, "low")  ~ "Low",
      TRUE ~ NA_character_
    ),
    Direction = dplyr::case_when(
      stringr::str_ends(taskCondition, "Acc") ~ "Faster",
      stringr::str_ends(taskCondition, "Dec") ~ "Slower",
      TRUE ~ NA_character_
    ),
    run = dplyr::if_else(ses == "ses1", 1L, 2L),
    # Condition also NA for ses2 Visual-first — same fix as long
    Condition = dplyr::case_when(
      !is.na(Condition) ~ Condition,
      TRUE              ~ NA_character_
    )
  ) %>%
  dplyr::left_join(
    summary_data %>% dplyr::select(id, Group),
    by = "id"
  ) %>%
  dplyr::mutate(
    Condition = dplyr::case_when(
      !is.na(Condition)                       ~ Condition,
      Group == "breath" & ses == "ses1"       ~ "breath",
      Group == "breath" & ses == "ses2"       ~ "breath",
      Group == "visual" & ses == "ses1"       ~ "visual",
      Group == "visual" & ses == "ses2"       ~ "breath",
      TRUE ~ NA_character_
    )
  )

# ---- 2c. Threshold data: reshape wide → long for H1/H2 ----
# summary_data has columns thresh_ses{1,2}_{Low,High}_{Faster,Slower}
# Reshape to one row per participant × session × salience × direction
threshold_long <- summary_data %>%
  dplyr::select(id, Group,
                starts_with("thresh_")) %>%
  tidyr::pivot_longer(
    cols      = starts_with("thresh_"),
    names_to  = c("ses", "Salience", "Direction"),
    names_pattern = "thresh_(ses[12])_(Low|High)_(Faster|Slower)",
    values_to = "Threshold"
  ) %>%
  dplyr::mutate(
    run       = dplyr::if_else(ses == "ses1", 1L, 2L),
    Salience  = factor(Salience,  levels = c("High", "Low")),
    Direction = factor(Direction, levels = c("Faster", "Slower"))
  ) %>%
  dplyr::filter(!is.na(Threshold))

cat(sprintf("Threshold long format: %d rows, %d participants\n",
            nrow(threshold_long), dplyr::n_distinct(threshold_long$id)))

# ---- 2d. Participant-level summary for H3, H6 ----
# Thresholds are condition-labelled using Group × ses mapping:
#   Breath-first: ses1 = breath condition, ses2 = breath condition
#   Visual-first: ses1 = visual condition, ses2 = breath condition
# This produces mean_breath_threshold and mean_visual_threshold
# that are comparable across groups, and overall_mean_threshold
# that averages all available data for H3.

threshold_by_condition <- threshold_long %>%
  dplyr::mutate(
    task_condition = dplyr::case_when(
      Group == "breath" & ses == "ses1" ~ "breath",
      Group == "breath" & ses == "ses2" ~ "breath",
      Group == "visual" & ses == "ses1" ~ "visual",
      Group == "visual" & ses == "ses2" ~ "breath",
      TRUE ~ NA_character_
    )
  )

participant_data <- summary_data %>%
  dplyr::left_join(
    # Mean threshold across all sessions and conditions (for H3)
    threshold_by_condition %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(
        overall_mean_threshold = mean(Threshold, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "id"
  ) %>%
  dplyr::left_join(
    # Condition-specific threshold means
    threshold_by_condition %>%
      dplyr::filter(!is.na(task_condition)) %>%
      dplyr::group_by(id, task_condition) %>%
      dplyr::summarise(
        mean_threshold = mean(Threshold, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      tidyr::pivot_wider(
        names_from  = task_condition,
        values_from = mean_threshold,
        names_prefix = "mean_thresh_"
      ),
    by = "id"
  ) %>%
  dplyr::mutate(
    # Mean confidence across all available sessions
    mean_confidence = rowMeans(
      dplyr::select(., mean_Confidence_ses1, mean_Confidence_ses2),
      na.rm = TRUE
    )
  )

cat(sprintf(
  "Participant summary: %d participants\n  overall threshold n=%d  breath n=%d  visual n=%d\n",
  nrow(participant_data),
  sum(!is.na(participant_data$overall_mean_threshold)),
  sum(!is.na(participant_data$mean_thresh_breath)),
  sum(!is.na(participant_data$mean_thresh_visual))
))
# Confirm the condition-specific columns exist
stopifnot(
  "mean_thresh_breath not found — check threshold_by_condition task_condition values" =
    "mean_thresh_breath" %in% names(participant_data),
  "mean_thresh_visual not found — check threshold_by_condition task_condition values" =
    "mean_thresh_visual" %in% names(participant_data)
)

# ---- 2e. Paced trial data for H4 (both conditions) ----
# Keep both breath and visual condition paced trials.
# H4C tests whether the change–arousal relationship differs
# between conditions, so both are needed in the join.
# Condition filtering happens inside the H4 models where relevant.
long_paced <- long %>%
  dplyr::filter(
    tolower(taskCondition) %in% c("highsalienceacc","highsaliencedec",
                                   "lowsalienceacc","lowsaliencedec",
                                   "breath"),
    !is.na(Arousal),
    !is.na(Change)
  )

cat(sprintf("Paced trials in long (both conditions): %d rows, %d participants\n",
            nrow(long_paced), dplyr::n_distinct(long_paced$id)))
cat(sprintf("  breath: %d  |  visual: %d\n",
            sum(long_paced$Condition == "breath", na.rm = TRUE),
            sum(long_paced$Condition == "visual", na.rm = TRUE)))


# =============================================================
#  3. BUILD PHYSIO REGIME SUBSETS
#  (applied to belt QC data for H4 sensitivity analysis)
# =============================================================

# ── Regime 2b exclusion: Breath-session poor synchrony ──────
qcSummary <- qcSummary %>%
  dplyr::mutate(
    breath_run = dplyr::if_else(first_condition == "Breath", 1L, 2L)
  )

qcSummary_breath <- qcSummary %>%
  dplyr::filter(run == breath_run)

regime2b_exclude <- qcSummary_breath %>%
  dplyr::filter(
    r_delta_change_all < 0.30,
    pct_correct_breaths < 50
  ) %>%
  dplyr::pull(id)

# ── Regime 3 exclusion: unusable belt ───────────────────────
regime3_exclude <- qcSummary %>%
  dplyr::filter(belt_quality == "unusable") %>%
  dplyr::pull(id) %>%
  unique()

message(sprintf("Regime 2b: %d participants excluded: %s",
                length(regime2b_exclude),
                paste(sort(regime2b_exclude), collapse = ", ")))
message(sprintf("Regime 3:  %d additional unusable-belt exclusions",
                length(setdiff(regime3_exclude, regime2b_exclude))))

# ── Trial-level QC flags (paced trials only) ────────────────
paced_conds <- c("breath","lowsalienceacc","lowsaliencedec",
                 "highsalienceacc","highsaliencedec")

qc_paced <- qcFull %>%
  dplyr::filter(
    trial_available == TRUE,
    tolower(condition) %in% paced_conds
  ) %>%
  dplyr::mutate(
    # Reconstruct task condition directly from qcFull fields.
    # Breath-first (first_condition=="Breath"): both runs = breath.
    # Visual-first (first_condition=="Visual"): run 1 = visual, run 2 = breath.
    Condition = dplyr::case_when(
      first_condition == "Breath" & run == 1 ~ "breath",
      first_condition == "Breath" & run == 2 ~ "breath",
      first_condition == "Visual" & run == 1 ~ "visual",
      first_condition == "Visual" & run == 2 ~ "breath",
      TRUE ~ NA_character_
    ),
    flag_cb  = correct_breaths == FALSE,
    flag_lag = lag_flag == TRUE,
    flag_dur = is.na(dur_b1) | dur_b1 < 2.0 | dur_b1 > 7.0,
    n_flags  = as.integer(flag_cb) + as.integer(flag_lag) +
               as.integer(flag_dur),
    regime2_bad = n_flags >= 2
  )

cat("\n--- Condition check in qc_paced ---\n")
print(qc_paced %>%
  dplyr::count(first_condition, run, Condition) %>%
  as.data.frame())

# ── Four regime datasets ─────────────────────────────────────
dat_r1  <- qc_paced
dat_r2  <- dplyr::filter(qc_paced, !regime2_bad)
dat_r2b <- dplyr::filter(qc_paced, !regime2_bad, !id %in% regime2b_exclude)
dat_r3  <- dplyr::filter(qc_paced, !regime2_bad,
                          !id %in% regime2b_exclude,
                          !id %in% regime3_exclude)

regime_datasets <- list(R1 = dat_r1, R2 = dat_r2,
                        R2b = dat_r2b, R3 = dat_r3)

cat("\n--- Regime sample sizes (paced trials) ---\n")
for (rname in names(regime_datasets)) {
  d <- regime_datasets[[rname]]
  cat(sprintf("  %s: %d trials, %d participants\n",
              rname, nrow(d), dplyr::n_distinct(d$id)))
}

# ── Join arousal from long file into each regime ─────────────
# trial_num in qcFull is sequential 1–80 across both runs.
# Trial in long resets within each taskCondition.
# Fix: both sides get row_number() within id × run — no reordering,
# just label each row's position as-is and join on that.
long_paced <- long_paced %>%
  dplyr::group_by(id, run) %>%
  dplyr::mutate(trial_within_run = dplyr::row_number()) %>%
  dplyr::ungroup()

add_arousal_keyed <- function(dat) {
  dat %>%
    dplyr::group_by(id, run) %>%
    dplyr::mutate(trial_within_run = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
      long_paced %>%
        dplyr::select(id, run, trial_within_run,
                      Arousal, Accuracy, Confidence, Change),
      by = c("id", "run", "trial_within_run")
    )
}

regime_arousal <- lapply(regime_datasets, add_arousal_keyed)

# Spot-check join quality
n_matched <- sum(!is.na(regime_arousal$R1$Arousal))
n_total   <- nrow(regime_arousal$R1)
cat(sprintf("\nJoin check (R1): %d / %d trials have Arousal (%.0f%%)\n",
            n_matched, n_total, 100 * n_matched / n_total))


# =============================================================
#  4. BEHAVIOURAL HYPOTHESES (H1, H2, H3, H5, H6)
#  Run once — physio regime does not apply
# =============================================================
cat("\n\n========================================\n")
cat("BEHAVIOURAL HYPOTHESES\n")
cat("========================================\n")

# ---- H1: Direction effect on threshold ----
# Prediction: Faster threshold < Slower threshold
# (acceleration easier to detect than deceleration)
# Model uses threshold_long with Direction as fixed effect
cat("\n--- H1: Threshold ~ Direction ---\n")
model_H1 <- lmer(
  Threshold ~ Direction + (1 | id),
  data = threshold_long,
  REML = FALSE
)
print(summary(model_H1)$coefficients)
cat(sprintf("  Direction effect (Slower vs Faster): β = %.3f, p = %.4f\n",
            fixef(model_H1)["DirectionSlower"],
            summary(model_H1)$coefficients["DirectionSlower", "Pr(>|t|)"]))

# ---- H2: Salience effect on threshold ----
# Prediction: High salience threshold < Low salience threshold
cat("\n--- H2: Threshold ~ Salience ---\n")
model_H2 <- lmer(
  Threshold ~ Salience + (1 | id),
  data = threshold_long,
  REML = FALSE
)
print(summary(model_H2)$coefficients)
cat(sprintf("  Salience effect (Low vs High): β = %.3f, p = %.4f\n",
            fixef(model_H2)["SalienceLow"],
            summary(model_H2)$coefficients["SalienceLow", "Pr(>|t|)"]))

# ---- H3: BCAT × MAIA correlations ----
# H3A: Confidence ~ MAIA_total (expect positive r)
# H3B: Threshold  ~ MAIA_total (expect null; BF01 ≥ 3)
# Both use all available data (both sessions averaged).
cat("\n--- H3A: Mean confidence (both sessions) ~ MAIA_total ---\n")
h3_data <- participant_data %>%
  dplyr::filter(!is.na(MAIA_total), !is.na(mean_confidence))

model_H3A <- cor.test(h3_data$mean_confidence, h3_data$MAIA_total)
print(model_H3A)

cat("\n--- H3B: Overall mean threshold (both sessions) ~ MAIA_total ---\n")
h3b_data <- participant_data %>%
  dplyr::filter(!is.na(MAIA_total), !is.na(overall_mean_threshold))

model_H3B_freq <- cor.test(h3b_data$overall_mean_threshold,
                            h3b_data$MAIA_total)
print(model_H3B_freq)

# Bayesian: expect BF01 ≥ 3 (evidence for null)
bf_H3B  <- BayesFactor::correlationBF(
  h3b_data$overall_mean_threshold,
  h3b_data$MAIA_total
)
# @bayesFactor$bf stores the log BF10; exp() converts to BF10.
# Error = 0% is correct — correlationBF uses an analytic solution
# (Jeffreys-beta* prior), not MCMC sampling, so there is no
# Monte Carlo error. Non-zero error only appears for sampled BFs
# (e.g. ttestBF, lmBF).
bf10 <- exp(bf_H3B@bayesFactor$bf)
bf01 <- 1 / bf10
cat(sprintf("  BF10 = %.3f  |  BF01 = %.3f\n", bf10, bf01))
cat(sprintf("  Criterion BF01 ≥ 3: %s\n",
            ifelse(bf01 >= 3, "MET — supports null", "NOT MET")))

# ---- H3 Exploratory: condition- and session-specific correlations ~ MAIA ----
# Four cells: ses1 × Breath-first, ses1 × Visual-first,
#             ses2 × Breath-first, ses2 × Visual-first
# For each cell: frequentist r + BF10 + BF01 for both threshold and confidence.
# BF01 ≥ 3 = evidence for null; BF10 ≥ 3 = evidence for alternative.
# Error = 0% in all cases (correlationBF uses analytic solution).
cat("\n--- H3 Exploratory: threshold and confidence ~ MAIA by ses × group ---\n")
cat(sprintf("  %-28s  %6s  %8s  %8s  %8s  %8s\n",
            "Cell", "n", "r", "p", "BF10", "BF01"))
cat(paste(rep("-", 78), collapse = ""), "\n")

# Helper: run frequentist + Bayesian correlation, print one row
report_cor <- function(label, x, y) {
  ok  <- !is.na(x) & !is.na(y)
  n   <- sum(ok)
  if (n < 5) {
    cat(sprintf("  %-28s  %6d  %8s  %8s  %8s  %8s\n",
                label, n, "—", "—", "—", "— (n too small)"))
    return(invisible(NULL))
  }
  fr  <- cor.test(x[ok], y[ok])
  bf  <- BayesFactor::correlationBF(x[ok], y[ok])
  bf10 <- exp(bf@bayesFactor$bf)
  bf01 <- 1 / bf10
  cat(sprintf("  %-28s  %6d  %+7.3f  %8.4f  %8.3f  %8.3f%s\n",
              label, n,
              fr$estimate, fr$p.value,
              bf10, bf01,
              ifelse(bf01 >= 3, " *H0", ifelse(bf10 >= 3, " *H1", ""))))
  invisible(list(r = fr$estimate, p = fr$p.value, bf10 = bf10, bf01 = bf01))
}

# Build a long table with ses × group × threshold and confidence
# Rejoin threshold_long (has per-ses values) with Group and MAIA
thresh_maia <- threshold_by_condition %>%
  dplyr::group_by(id, ses, task_condition) %>%
  dplyr::summarise(mean_thresh = mean(Threshold, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::left_join(
    summary_data %>% dplyr::select(id, Group, MAIA_total,
                                    mean_Confidence_ses1,
                                    mean_Confidence_ses2),
    by = "id"
  ) %>%
  dplyr::mutate(
    mean_conf = dplyr::if_else(ses == "ses1",
                               mean_Confidence_ses1,
                               mean_Confidence_ses2)
  ) %>%
  dplyr::filter(!is.na(MAIA_total))

cat("\n  THRESHOLD ~ MAIA\n")
for (tc in c("breath", "visual")) {
  for (g in c("breath", "visual")) {
    sub <- thresh_maia %>%
      dplyr::filter(task_condition == tc, Group == g, !is.na(mean_thresh))
    label <- sprintf("%s cond | %s-first", tc, g)
    report_cor(label, sub$mean_thresh, sub$MAIA_total)
  }
}

cat("\n  CONFIDENCE ~ MAIA\n")
for (tc in c("breath", "visual")) {
  for (g in c("breath", "visual")) {
    sub <- thresh_maia %>%
      dplyr::filter(task_condition == tc, Group == g, !is.na(mean_conf))
    label <- sprintf("%s cond | %s-first", tc, g)
    report_cor(label, sub$mean_conf, sub$MAIA_total)
  }
}

# Also: condition-level means (breath vs visual, collapsing across ses)
cat("\n  CONDITION-LEVEL (both sessions pooled)\n")
for (cond_label in c("breath", "visual")) {
  col <- paste0("mean_thresh_", cond_label)
  sub <- participant_data %>% dplyr::filter(!is.na(MAIA_total), !is.na(.data[[col]]))
  report_cor(sprintf("threshold | %s cond", cond_label),
             sub[[col]], sub$MAIA_total)
}
report_cor("confidence | overall",
           participant_data %>% dplyr::filter(!is.na(MAIA_total)) %>% dplyr::pull(mean_confidence),
           participant_data %>% dplyr::filter(!is.na(MAIA_total)) %>% dplyr::pull(MAIA_total))

# Breath vs visual threshold correlation
cat("\n  CROSS-CONDITION\n")
trt_cond <- participant_data %>%
  dplyr::filter(!is.na(mean_thresh_breath), !is.na(mean_thresh_visual))
r_bv <- cor.test(trt_cond$mean_thresh_breath, trt_cond$mean_thresh_visual)
cat(sprintf("  Breath vs visual threshold (n=%d): r = %.3f  p = %.4f\n",
            nrow(trt_cond), r_bv$estimate, r_bv$p.value))

# ---- H5: Test trial salience effect on accuracy ----
# Prediction: High salience → higher detection accuracy on test block
# Uses test-phase data only; Breath-condition sessions only
cat("\n--- H5: Accuracy ~ Salience (test trials, Breath condition) ---\n")
test_breath <- test %>%
  dplyr::filter(Condition == "breath", !is.na(Salience))

model_H5 <- lmer(
  Accuracy ~ Salience + (1 | id),
  data = test_breath,
  REML = FALSE
)
print(summary(model_H5)$coefficients)
cat(sprintf("  Salience effect (Low vs High): β = %.3f, p = %.4f\n",
            fixef(model_H5)["SalienceLow"],
            summary(model_H5)$coefficients["SalienceLow", "Pr(>|t|)"]))

# ---- H6: Test-retest reliability ----
# Prediction: r ≥ .70 overall; lower in Visual-first than Breath-first.
# Test-retest is computed between sessions, which means between
# breath and visual condition thresholds (each participant contributes
# one of each). The correlation between breath and visual thresholds
# also quantifies cross-condition transfer.
cat("\n--- H6: Test-retest reliability ---\n")

# Overall: ses1 vs ses2 regardless of condition
trt_data <- summary_data %>%
  dplyr::left_join(
    threshold_long %>%
      dplyr::group_by(id, ses) %>%
      dplyr::summarise(mean_thresh = mean(Threshold, na.rm = TRUE),
                       .groups = "drop") %>%
      tidyr::pivot_wider(names_from = ses, values_from = mean_thresh,
                         names_prefix = "thresh_"),
    by = "id"
  ) %>%
  dplyr::filter(!is.na(thresh_ses1), !is.na(thresh_ses2))

trt_overall <- cor.test(trt_data$thresh_ses1, trt_data$thresh_ses2)
cat(sprintf("  Overall ses1 vs ses2: r = %.3f [%.3f, %.3f]  p = %.4f  n = %d\n",
            trt_overall$estimate,
            trt_overall$conf.int[1], trt_overall$conf.int[2],
            trt_overall$p.value, nrow(trt_data)))
cat(sprintf("  Target r ≥ .70: %s\n",
            ifelse(trt_overall$estimate >= .70, "MET", "NOT MET")))

# By group — tests whether order effect moderates reliability
for (grp in c("breath", "visual")) {
  sub <- trt_data %>% dplyr::filter(Group == grp)
  r   <- cor.test(sub$thresh_ses1, sub$thresh_ses2)
  cat(sprintf("  %s-first: r = %.3f [%.3f, %.3f]  n = %d\n",
              grp, r$estimate, r$conf.int[1], r$conf.int[2], nrow(sub)))
}

# Exploratory: breath vs visual threshold correlation
# (cross-condition reliability — are participants who have low breath
# thresholds also better at detecting visual changes?)
cat("\n  Exploratory — breath vs visual threshold correlation:\n")
bv_data <- participant_data %>%
  dplyr::filter(!is.na(mean_thresh_breath), !is.na(mean_thresh_visual))
r_bv2 <- cor.test(bv_data$mean_thresh_breath, bv_data$mean_thresh_visual)
cat(sprintf("  r = %.3f [%.3f, %.3f]  p = %.4f  n = %d\n",
            r_bv2$estimate, r_bv2$conf.int[1], r_bv2$conf.int[2],
            r_bv2$p.value, nrow(bv_data)))


# =============================================================
#  5. BELT-INFORMED HYPOTHESES (H4A, H4B, H4C) × 4 REGIMES
# =============================================================
cat("\n\n========================================\n")
cat("H4: AROUSAL × RESPIRATORY CHANGE\n")
cat("========================================\n")

# ── Helper: run all H4 models on one regime ──────────────────
run_H4_models <- function(dat, regime_label) {

  dat_h4 <- dat %>% dplyr::filter(!is.na(Arousal), !is.na(Change))

  # H4A: Arousal ~ Change (breath condition only)
  # Negative β expected (faster = negative Change → higher arousal).
  # Belt-informed: uses qc_paced × long join. Breath condition only
  # since qcFull only processes Breath-condition sessions.
  m_H4A <- lmer(Arousal ~ Change + (1 | id),
                data = dat_h4,
                REML = FALSE)

  # H4B: Arousal ~ Change × Accuracy (breath condition only)
  # Preregistered moderation: change–arousal link stronger when
  # participant correctly detected the change (Accuracy = 1).
  m_H4B <- lmer(Arousal ~ Change * Accuracy + (1 | id),
                data = dat_h4 %>% dplyr::filter(!is.na(Accuracy)),
                REML = FALSE)

  # H4C is run separately on long (behavioural) — see below.
  # Belt QC pipeline only processed Breath-condition sessions, so
  # Visual-condition trials do not exist in qc_paced / regime_arousal.
  # H4C is therefore purely behavioural and does not vary by regime.

  extract_fe <- function(model, hyp) {
    broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
      dplyr::mutate(hypothesis = hyp, regime = regime_label)
  }

  dplyr::bind_rows(
    extract_fe(m_H4A, "H4A"),
    extract_fe(m_H4B, "H4B")
  )
}

# ── Run across all four regimes ───────────────────────────────
h4_results <- purrr::map_dfr(
  names(regime_arousal),
  ~ run_H4_models(regime_arousal[[.x]], .x)
)

# ── Primary results (Regime 1) ────────────────────────────────
cat("\n--- H4 primary results (Regime 1) ---\n")
print(
  h4_results %>%
    dplyr::filter(regime == "R1") %>%
    dplyr::select(hypothesis, term, estimate, std.error,
                  statistic, p.value, conf.low, conf.high) %>%
    as.data.frame(),
  digits = 3, row.names = FALSE
)

# ── H4C: behavioural model (both conditions, no regime comparison) ────────
# Run on long directly. qcFull only contains Breath-condition sessions,
# so H4C cannot be tested via regime_arousal. Change here is the
# staircase-log intended change, not belt-measured.
cat("\n--- H4C: Arousal ~ Change × Condition (behavioural, all participants) ---\n")
long_h4c <- long %>%
  dplyr::filter(
    tolower(taskCondition) %in% c("highsalienceacc","highsaliencedec",
                                   "lowsalienceacc","lowsaliencedec"),
    !is.na(Arousal), !is.na(Change), !is.na(Condition)
  )
cat(sprintf("  n = %d trials, %d participants; Condition: %s\n",
            nrow(long_h4c), dplyr::n_distinct(long_h4c$id),
            paste(sort(unique(long_h4c$Condition)), collapse = ", ")))

model_H4C <- lmer(Arousal ~ Change * Condition + (1 | id),
                  data = long_h4c, REML = FALSE)
print(broom.mixed::tidy(model_H4C, effects = "fixed", conf.int = TRUE),
      n = Inf)
cat(sprintf("  H4C interaction (visual vs breath): β = %.3f, p = %.4f\n",
            fixef(model_H4C)["Change:Conditionvisual"],
            summary(model_H4C)$coefficients["Change:Conditionvisual","Pr(>|t|)"]))


# =============================================================
#  6. BENJAMINI-HOCHBERG CORRECTION
# =============================================================
cat("\n\n========================================\n")
cat("BENJAMINI-HOCHBERG CORRECTION\n")
cat("========================================\n")

# ── Family 1: Primary main effects ───────────────────────────
h1_p <- summary(model_H1)$coefficients["DirectionSlower",  "Pr(>|t|)"]
h2_p <- summary(model_H2)$coefficients["SalienceLow",      "Pr(>|t|)"]
h4a_p <- h4_results %>%
  dplyr::filter(regime == "R1", hypothesis == "H4A", term == "Change") %>%
  dplyr::pull(p.value)
h4c_int_p <- summary(model_H4C)$coefficients[
  "Change:Conditionvisual", "Pr(>|t|)"]
h5_p <- summary(model_H5)$coefficients["SalienceLow", "Pr(>|t|)"]

family1 <- tibble::tibble(
  hypothesis = c("H1", "H2", "H4A", "H4C_interaction", "H5"),
  term       = c("DirectionSlower", "SalienceLow", "Change",
                 "Change:Conditionvisual", "SalienceLow"),
  p_raw      = c(h1_p, h2_p, h4a_p, h4c_int_p, h5_p)
) %>%
  dplyr::filter(!is.na(p_raw)) %>%
  dplyr::mutate(p_BH = p.adjust(p_raw, method = "BH"),
                sig_BH = p_BH < .05)

cat("\nFamily 1 — main effects:\n")
print(as.data.frame(family1), digits = 4, row.names = FALSE)

# ── Family 2: Interaction terms (H4B only — H4C is behavioural standalone) ─
family2 <- h4_results %>%
  dplyr::filter(
    regime == "R1",
    hypothesis == "H4B",
    stringr::str_detect(term, ":")
  ) %>%
  dplyr::select(hypothesis, term, p_raw = p.value) %>%
  dplyr::bind_rows(
    tibble::tibble(
      hypothesis = "H4C",
      term       = "Change:Conditionvisual",
      p_raw      = h4c_int_p
    )
  ) %>%
  dplyr::mutate(p_BH = p.adjust(p_raw, method = "BH"),
                sig_BH = p_BH < .05)

cat("\nFamily 2 — interaction terms:\n")
print(as.data.frame(family2), digits = 4, row.names = FALSE)


# =============================================================
#  7. REGIME COMPARISON TABLE (H4)
# =============================================================
cat("\n\n========================================\n")
cat("REGIME COMPARISON\n")
cat("========================================\n")

key_terms <- c("Change", "Change:Accuracy")

regime_comparison <- h4_results %>%
  dplyr::filter(term %in% key_terms) %>%
  dplyr::select(hypothesis, term, regime,
                estimate, std.error, statistic,
                p.value, conf.low, conf.high) %>%
  dplyr::arrange(hypothesis, term, regime)

cat("\nEffect estimates across regimes:\n")
print(as.data.frame(regime_comparison), digits = 3, row.names = FALSE)

# Agreement summary
regime_agreement <- regime_comparison %>%
  dplyr::group_by(hypothesis, term) %>%
  dplyr::summarise(
    consistent_direction = dplyr::n_distinct(sign(estimate)) == 1,
    n_sig_raw            = sum(p.value < .05, na.rm = TRUE),
    estimate_range       = sprintf("[%.3f, %.3f]",
                                   min(estimate), max(estimate)),
    .groups = "drop"
  )
cat("\nRegime agreement:\n")
print(as.data.frame(regime_agreement), row.names = FALSE)


# =============================================================
#  8. EXPLORATORY BELT ANALYSES (clearly labelled)
# =============================================================
cat("\n\n========================================\n")
cat("EXPLORATORY BELT ANALYSES\n")
cat("(not preregistered — interpret cautiously)\n")
cat("========================================\n")

# ---- 8a. Does observed belt change predict arousal beyond intended? ----
cat("\n--- Exploratory 1: Arousal ~ Change + dur_b1_vs_b4 ---\n")
explore_belt <- regime_arousal$R2 %>%
  dplyr::filter(!is.na(Arousal), !is.na(dur_b1_vs_b4), !is.na(Change))

m_ex1 <- lmer(Arousal ~ Change + dur_b1_vs_b4 + (1 | id),
              data = explore_belt, REML = FALSE)
print(broom.mixed::tidy(m_ex1, effects = "fixed", conf.int = TRUE),
      n = Inf)

# ---- 8b. Belt quality as moderator ----
cat("\n--- Exploratory 2: Arousal ~ Change × belt_quality ---\n")
m_ex2 <- lmer(Arousal ~ Change * belt_quality + (1 | id),
              data = regime_arousal$R2 %>% dplyr::filter(!is.na(Arousal)),
              REML = FALSE)
print(broom.mixed::tidy(m_ex2, effects = "fixed", conf.int = TRUE),
      n = Inf)

# ---- 8c. MAIA moderation ----
cat("\n--- Exploratory 3: Arousal ~ Change × MAIA_total ---\n")
explore_maia <- regime_arousal$R2 %>%
  dplyr::left_join(
    summary_data %>% dplyr::select(id, MAIA_total),
    by = "id"
  ) %>%
  dplyr::filter(!is.na(Arousal), !is.na(MAIA_total)) %>%
  dplyr::mutate(maia_c = scale(MAIA_total, scale = FALSE)[, 1])

m_ex3 <- lmer(Arousal ~ Change * maia_c + (1 | id),
              data = explore_maia, REML = FALSE)
print(broom.mixed::tidy(m_ex3, effects = "fixed", conf.int = TRUE),
      n = Inf)

# ---- 8d. Decomposing intended change, observed compliance, and detection ----
# Compares four nested models to estimate relative contribution of:
#   (1) Intended change (pacer)
#   (2) Observed physiological compliance (belt)
#   (3) Conscious detection (Accuracy)
# Uses Regime 2 (trial-level QC applied, belt signal reliable).
# Standardised predictors allow direct β comparison across models.
cat("\n--- Exploratory 4: Decomposing intended change vs compliance vs detection ---\n")

explore_decomp <- regime_arousal$R2 %>%
  dplyr::filter(!is.na(Arousal), !is.na(Change),
                !is.na(dur_b1_vs_b4), !is.na(Accuracy)) %>%
  dplyr::mutate(
    Change_z      = scale(Change)[, 1],
    dur_b1_vs_b4_z = scale(dur_b1_vs_b4)[, 1]
  )

# M1: intended change only (baseline)
m_d1 <- lmer(Arousal ~ Change_z + (1 | id),
             data = explore_decomp, REML = FALSE)

# M2: observed belt change only
m_d2 <- lmer(Arousal ~ dur_b1_vs_b4_z + (1 | id),
             data = explore_decomp, REML = FALSE)

# M3: both change sources (is belt change unique beyond intended?)
m_d3 <- lmer(Arousal ~ Change_z + dur_b1_vs_b4_z + (1 | id),
             data = explore_decomp, REML = FALSE)

# M4: add conscious detection
m_d4 <- lmer(Arousal ~ Change_z + dur_b1_vs_b4_z + Accuracy + (1 | id),
             data = explore_decomp, REML = FALSE)

# M5: does detection interact with belt compliance?
# Tests whether the belt→arousal pathway is amplified by detection
m_d5 <- lmer(Arousal ~ Change_z + dur_b1_vs_b4_z * Accuracy + (1 | id),
             data = explore_decomp, REML = FALSE)

tab_model(m_d1, m_d2, m_d3, m_d4, m_d5)
anova(m_d1, m_d2, m_d3, m_d4, m_d5)

# Print model comparison table
cat("\n  Model comparison (AIC, marginal R²):\n")
decomp_models <- list(
  "M1: Intended change"              = m_d1,
  "M2: Observed belt change"         = m_d2,
  "M3: Both change sources"          = m_d3,
  "M4: + Conscious detection"        = m_d4,
  "M5: + Belt × Detection interact." = m_d5
)

for (nm in names(decomp_models)) {
  m   <- decomp_models[[nm]]
  r2  <- MuMIn::r.squaredGLMM(m)
  cat(sprintf("  %-38s  AIC = %7.1f  R²m = %.3f  R²c = %.3f\n",
              nm, AIC(m), r2[1, "R2m"], r2[1, "R2c"]))
}

cat("\n  Fixed effects — M5 (full model):\n")
print(broom.mixed::tidy(m_d5, effects = "fixed", conf.int = TRUE), n = Inf)

# ---- 8e. Non-conscious compliance: do misses still breathe correctly? ----
# Key question: when participants MISS detection (Accuracy=0), do they still
# physically change their breathing in the correct direction? And does that
# physiological change still predict subjective arousal?
# This directly tests whether the belt→arousal pathway operates
# without conscious detection — the core misattribution mechanism.
cat("\n--- Exploratory 5: Non-conscious compliance on missed trials ---\n")

explore_miss <- regime_arousal$R2 %>%
  dplyr::filter(!is.na(Accuracy), !is.na(direction_correct),
                delta != 0)  # change trials only

# 5a. Direction correct rate by detection outcome
cat("\n  5a. Direction correct rate by detection accuracy:\n")
dir_by_acc <- explore_miss %>%
  dplyr::group_by(Accuracy) %>%
  dplyr::summarise(
    n              = n(),
    pct_dir_correct = mean(direction_correct, na.rm = TRUE) * 100,
    mean_belt_change = mean(dur_b1_vs_b4, na.rm = TRUE),
    sd_belt_change   = sd(dur_b1_vs_b4, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(dir_by_acc), digits = 3, row.names = FALSE)

# Is direction correct rate on misses above chance (50%)?
miss_dir <- explore_miss %>%
  dplyr::filter(Accuracy == 0, !is.na(direction_correct))
binom_test <- binom.test(
  x = sum(miss_dir$direction_correct, na.rm = TRUE),
  n = sum(!is.na(miss_dir$direction_correct)),
  p = 0.5
)
cat(sprintf("\n  Binomial test: direction correct on misses = %.1f%%\n",
            mean(miss_dir$direction_correct, na.rm=TRUE)*100))
cat(sprintf("  vs chance (50%%): p = %.4f  95%% CI [%.3f, %.3f]\n",
            binom_test$p.value,
            binom_test$conf.int[1], binom_test$conf.int[2]))

# 5b. Belt change magnitude on hits vs misses
# If similar magnitude, the physiological change happened regardless
cat("\n  5b. Belt change magnitude (dur_b1_vs_b4) by detection accuracy:\n")
belt_by_acc <- explore_miss %>%
  dplyr::filter(!is.na(dur_b1_vs_b4)) %>%
  dplyr::group_by(Accuracy) %>%
  dplyr::summarise(
    n      = n(),
    mean   = mean(dur_b1_vs_b4),
    sd     = sd(dur_b1_vs_b4),
    median = median(dur_b1_vs_b4),
    .groups = "drop"
  )
print(as.data.frame(belt_by_acc), digits = 3, row.names = FALSE)

# Wilcoxon test: is belt change magnitude different between hits and misses?
hit_belt  <- explore_miss %>% dplyr::filter(Accuracy == 1) %>%
  dplyr::pull(dur_b1_vs_b4) %>% na.omit()
miss_belt <- explore_miss %>% dplyr::filter(Accuracy == 0) %>%
  dplyr::pull(dur_b1_vs_b4) %>% na.omit()
wt_belt <- wilcox.test(hit_belt, miss_belt)
cat(sprintf("  Wilcoxon W = %.0f  p = %.4f\n", wt_belt$statistic, wt_belt$p.value))

# 5c. Does belt change predict arousal on MISSED trials specifically?
# The critical test: physiological arousal → subjective arousal
# without conscious detection as a mediator
cat("\n  5c. Belt change → Arousal separately for hits and misses:\n")
for (acc_val in c(0, 1)) {
  label <- if (acc_val == 0) "MISSES (Accuracy=0)" else "HITS   (Accuracy=1)"
  sub   <- regime_arousal$R2 %>%
    dplyr::filter(Accuracy == acc_val,
                  !is.na(Arousal), !is.na(dur_b1_vs_b4),
                  delta != 0)
  m <- lmer(Arousal ~ dur_b1_vs_b4 + (1 | id),
            data = sub, REML = FALSE)
  coef_row <- summary(m)$coefficients["dur_b1_vs_b4", ]
  cat(sprintf("  %s  n=%d  β = %+.3f  SE = %.3f  p = %.4f\n",
              label, nrow(sub),
              coef_row["Estimate"],
              coef_row["Std. Error"],
              coef_row["Pr(>|t|)"]))
}

# 5d. Full model: belt change predicts arousal controlling for detection
# and the interaction — does detection amplify the belt→arousal link?
cat("\n  5d. Arousal ~ dur_b1_vs_b4 × Accuracy (change trials only):\n")
m_noncon <- lmer(Arousal ~ dur_b1_vs_b4 * Accuracy + (1 | id),
                 data = regime_arousal$R2 %>%
                   dplyr::filter(!is.na(Arousal), !is.na(dur_b1_vs_b4),
                                 !is.na(Accuracy), delta != 0),
                 REML = FALSE)
print(broom.mixed::tidy(m_noncon, effects = "fixed", conf.int = TRUE),
      n = Inf)
cat("  Interpretation: dur_b1_vs_b4 main effect = belt→arousal WITHOUT detection\n")
cat("  Interaction = additional amplification WHEN detection occurs\n")


# =============================================================
#  9. FIGURES
# =============================================================

# ---- F1: Regime comparison — key H4 effects ----
fig_regime <- regime_comparison %>%
  dplyr::mutate(
    regime     = factor(regime, levels = c("R1","R2","R2b","R3")),
    hypothesis = factor(hypothesis),
    term_label = dplyr::recode(term,
      "Change"             = "Main: Change",
      "Change:Accuracy"    = "Interaction: Change × Accuracy",
      "Change:Conditionvisual" = "Interaction: Change × Group"
    )
  ) %>%
  ggplot(aes(x = regime, y = estimate,
             ymin = conf.low, ymax = conf.high,
             colour = regime)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.7) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  facet_wrap(~ term_label, scales = "free_y", nrow = 1) +
  scale_colour_viridis_d(option = "D", end = 0.85, guide = "none") +
  labs(
    title    = "H4: Effect stability across exclusion regimes",
    subtitle = "β ± 95% CI; dashed line = 0; Regime 1 = primary",
    x = "Exclusion regime",
    y = "β (effect on Arousal)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title  = element_text(face = "bold"),
        strip.text  = element_text(face = "bold", size = 9),
        panel.spacing = unit(1.2, "lines"))

ggsave(paste0(analysisPath, "S5_H4_regime_comparison.png"),
       fig_regime, width = 13, height = 5, dpi = 120)

# ---- F2: Arousal by Change — raw scatter with trend ----
fig_arousal <- regime_arousal$R1 %>%
  dplyr::filter(!is.na(Arousal), !is.na(Change)) %>%
  dplyr::left_join(summary_data %>% dplyr::select(id, Group), by = "id") %>%
  dplyr::mutate(
    Salience = dplyr::case_when(
      stringr::str_starts(tolower(condition), "high") ~ "High",
      stringr::str_starts(tolower(condition), "low")  ~ "Low",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Salience)) %>%
  ggplot(aes(x = Change, y = Arousal)) +
  geom_point(alpha = 0.06, size = 0.7, colour = "#2c7be5") +
  geom_smooth(method = "lm", colour = "#e63946",
              linewidth = 1.0, se = TRUE) +
  facet_grid(Salience ~ Group,
             labeller = labeller(
               Group = c(breath = "Breath-first",
                         visual = "Visual-first"))) +
  scale_x_continuous(
    breaks = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1\n(Faster)", "-0.5", "0\n(No change)",
               "0.5", "1\n(Slower)")
  ) +
  labs(
    title    = "H4A: Arousal by respiratory change (Regime 1)",
    subtitle = "Each point = one trial; line = OLS trend",
    x        = "Change (negative = faster, positive = slower)",
    y        = "Arousal rating (1–6)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

ggsave(paste0(analysisPath, "S5_H4A_arousal_change.png"),
       fig_arousal, width = 10, height = 7, dpi = 120)

# ---- F3: Threshold by Direction and Salience (H1, H2) ----
fig_threshold <- threshold_long %>%
  ggplot(aes(x = Direction, y = Threshold, fill = Salience)) +
  geom_violin(alpha = 0.35, position = position_dodge(0.8),
              linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.7,
               position = position_dodge(0.8)) +
  scale_fill_manual(values = c("High" = "#2dc653", "Low" = "#f4a261")) +
  labs(
    title = "H1/H2: Threshold by Direction and Salience",
    x = "Direction", y = "Threshold (level)",
    fill = "Salience"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(paste0(analysisPath, "S5_H1H2_thresholds.png"),
       fig_threshold, width = 8, height = 5, dpi = 120)

message("\nanalysis_study5.R complete.")
