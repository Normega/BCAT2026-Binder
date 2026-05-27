# ============================================================
# analysis_study5_exploratory.R
# Study 5 exploratory MAIA and H3B analyses.
# Sources after analysis_maia.R; uses s5s, s5_long_breath.
# Outputs: s5e_h3b_threshold_breakdown.csv,
#          s5e_maia_subscale_cors.csv, s5e_alexithymia_dissociation.csv,
#          s5e_group_awareness.csv, s5e_training_contrast.csv,
#          s5e_awareness_change.csv
# ============================================================

source(file.path(ANALYSIS_DIR, "theme_bcat.R"))

s5e_fig_dir <- file.path(FIG_DIR, "Study5_Exploratory")
dir.create(s5e_fig_dir, showWarnings = FALSE, recursive = TRUE)


# s5s and s5_long_breath loaded by MainAnalysis.R
summary_data <- s5s
# Build long-format threshold: one row per participant × session ×
# salience × direction. This enables mixed models over all conditions.
thresh_long <- summary_data |>
  dplyr::select(id, Group, MAIA_total_z,
                dplyr::matches("^thresh_ses[12]_"),
                dplyr::matches("^MAIA_"),
                mean_Confidence_ses1, mean_Confidence_ses2,
                Awareness_ses1, Awareness_ses2,
                MSES_selfdoubt) |>
  tidyr::pivot_longer(
    cols          = dplyr::matches("^thresh_ses[12]_"),
    names_to      = c("ses", "Salience", "Direction"),
    names_pattern = "thresh_(ses[12])_([^_]+)_(.*)",
    values_to     = "Threshold"
  ) |>
  dplyr::mutate(
    Salience  = factor(Salience,  levels = c("Low",    "High")),
    Direction = factor(Direction, levels = c("Faster", "Slower")),
    ses       = factor(ses),
    Condition = dplyr::if_else(
      (Group == "breath") |
        (Group == "visual" & ses == "ses2"),
      "breath", "visual"
    ),
    Condition = factor(Condition, levels = c("breath", "visual"))
  ) |>
  dplyr::filter(!is.na(Threshold))

message(sprintf("Threshold long: %d rows, %d participants",
                nrow(thresh_long), dplyr::n_distinct(thresh_long$id)))


# ============================================================
# E1. H3B THREE-LEVEL BREAKDOWN
#
# Each level tests MAIA_total_z ~ threshold in a progressively
# more granular way. All tests at this level are exploratory.
# The primary H3B test (overall mean threshold) lives in unified.
# ============================================================
cat("\n\n========================================\n")
cat("E1: H3B THRESHOLD ~ MAIA — LEVEL BREAKDOWN\n")
cat("(All exploratory — not preregistered)\n")
cat("========================================\n")

# -- LEVEL 1: Overall (reference — matches unified) ------------
overall_thresh <- summary_data |>
  dplyr::mutate(
    overall_mean_threshold = rowMeans(
      dplyr::across(dplyr::matches("^thresh_ses1_")), na.rm = TRUE)
  )
cat("\nLevel 1 — Overall mean threshold (matches unified H3B primary):\n")
e1_overall <- cor.test(overall_thresh$overall_mean_threshold,
                       overall_thresh$MAIA_total_z)
cat(sprintf("  %s  N=%d\n", fmt_cor(e1_overall), nrow(
  overall_thresh |> dplyr::filter(!is.na(overall_mean_threshold),
                                   !is.na(MAIA_total_z)))))

# BF₀₁ for null
e1_bf <- tryCatch(
  BayesFactor::correlationBF(
    overall_thresh$overall_mean_threshold,
    overall_thresh$MAIA_total,
    iterations = 10000
  ), error = function(e) NULL)
if (!is.null(e1_bf)) {
  bf10 <- BayesFactor::extractBF(e1_bf)$bf
  cat(sprintf("  BF01 = %.3f (%s)\n", 1/bf10,
              dplyr::case_when(
                1/bf10 >= 10 ~ "strong null evidence",
                1/bf10 >= 3  ~ "moderate null evidence",
                TRUE         ~ "anecdotal")))
}

# -- LEVEL 2: Simple effects (Salience and Direction margins) --
cat("\nLevel 2 — Simple effects:\n")

level2_results <- purrr::map_dfr(
  list(
    list(var="Salience",  val="High",   col_regex="^thresh_ses1_High_"),
    list(var="Salience",  val="Low",    col_regex="^thresh_ses1_Low_"),
    list(var="Direction", val="Faster", col_regex="^thresh_ses1_.*_Faster"),
    list(var="Direction", val="Slower", col_regex="^thresh_ses1_.*_Slower")
  ),
  function(spec) {
    d <- summary_data |>
      dplyr::mutate(
        thresh_mean = rowMeans(
          dplyr::across(dplyr::matches(spec$col_regex)), na.rm = TRUE)
      ) |>
      dplyr::filter(!is.na(thresh_mean), !is.na(MAIA_total_z))
    ct <- cor.test(d$thresh_mean, d$MAIA_total_z)
    tibble::tibble(
      factor    = spec$var,
      level     = spec$val,
      r         = ct$estimate,
      ci_lower  = ct$conf.int[1],
      ci_upper  = ct$conf.int[2],
      p_value   = ct$p.value,
      n         = nrow(d)
    )
  }
)
print(as.data.frame(level2_results), digits = 3, row.names = FALSE)

# -- LEVEL 3: All 4 Salience × Direction cells -----------------
cat("\nLevel 3 — Salience × Direction cells:\n")

thresh_cols <- c("thresh_ses1_High_Faster", "thresh_ses1_High_Slower",
                 "thresh_ses1_Low_Faster",  "thresh_ses1_Low_Slower")

level3_results <- purrr::map_dfr(thresh_cols, function(col) {
  d <- summary_data |>
    dplyr::filter(!is.na(.data[[col]]), !is.na(MAIA_total_z))
  ct <- cor.test(d[[col]], d$MAIA_total_z)
  tibble::tibble(
    condition = col,
    r         = ct$estimate,
    ci_lower  = ct$conf.int[1],
    ci_upper  = ct$conf.int[2],
    p_value   = ct$p.value,
    n         = nrow(d)
  )
})
print(as.data.frame(level3_results), digits = 3, row.names = FALSE)

# -- LEVEL 4: Condition and session breakdown ------------------
cat("\nLevel 4 — Condition (breath/visual) and session:\n")

level4_results <- thresh_long |>
  dplyr::filter(!is.na(MAIA_total_z), !is.na(Threshold)) |>
  dplyr::group_by(id, Condition, ses, MAIA_total_z) |>
  dplyr::summarise(mean_threshold = mean(Threshold, na.rm = TRUE), .groups = "drop") |>
  dplyr::group_by(Condition, ses) |>
  dplyr::summarise(
    r       = cor(mean_threshold, MAIA_total_z, use = "complete.obs"),
    n_parts = dplyr::n_distinct(id),
    .groups = "drop"
  )
print(as.data.frame(level4_results), digits = 3, row.names = FALSE)

# Mixed model: MAIA_total_z as predictor of threshold across all
# conditions simultaneously, with random participant intercept
m_e1_full <- lme4::lmer(
  Threshold ~ MAIA_total_z * Salience * Direction +
    Condition + ses + (1 | id),
  data = thresh_long |> dplyr::filter(!is.na(MAIA_total_z)),
  REML = FALSE
)

m_e1_null <- lme4::lmer(
  Threshold ~ Salience * Direction + Condition + ses + (1 | id),
  data = thresh_long |> dplyr::filter(!is.na(MAIA_total_z)),
  REML = FALSE
)

cat("\nLRT: MAIA_total_z in threshold model across all conditions:\n")
print(anova(m_e1_null, m_e1_full))

# Save level breakdown
e1_all <- dplyr::bind_rows(
  level2_results |>
    dplyr::mutate(level = "Simple effects",
                  condition = paste(factor, level)),
  level3_results |>
    dplyr::mutate(level = "Cells", factor = "Salience×Direction")
)
readr::write_csv(e1_all,
                 file.path(RESULTS_DIR, "s5e_h3b_threshold_breakdown.csv"))
message("Saved: s5e_h3b_threshold_breakdown.csv")


# ============================================================
# E2. MAIA SUBSCALE CORRELATIONS
#
# Does any particular subscale correlate with threshold or
# confidence? Tests the MAIA-as-metacognition hypothesis at
# a finer grain: sensibility facets should correlate with
# confidence (Noticing, Attention-Regulation) while the
# non-distraction facets (whose reversal we just fixed)
# are hypothesized to be less predictive.
# ============================================================
cat("\n\n========================================\n")
cat("E2: MAIA SUBSCALE CORRELATIONS\n")
cat("========================================\n")

subscales <- c("MAIA_Noticing", "MAIA_NotDistracting", "MAIA_NotWorrying",
               "MAIA_AttentionReg", "MAIA_EmoAware", "MAIA_SelfReg",
               "MAIA_BodyListen", "MAIA_Trusting")

# Z-scored versions (from standardise_maia)
subscale_z_cols <- paste0(subscales, "_z")
available_z <- intersect(subscale_z_cols, names(summary_data))

outcomes_sum <- list(
  threshold  = overall_thresh |>
    dplyr::select(id, value = overall_mean_threshold),
  confidence = summary_data |>
    dplyr::mutate(value = mean_Confidence_ses1) |>
    dplyr::select(id, value)
)

subscale_cors <- purrr::map_dfr(names(outcomes_sum), function(outcome) {
  d_out <- outcomes_sum[[outcome]]
  purrr::map_dfr(available_z, function(sc) {
    sc_raw <- sub("_z$", "", sc)
    d <- summary_data |>
      dplyr::select(id, maia = dplyr::all_of(sc)) |>
      dplyr::left_join(d_out, by = "id") |>
      tidyr::drop_na()
    if (nrow(d) < 10) return(NULL)
    ct <- cor.test(d$value, d$maia)
    tibble::tibble(
      outcome   = outcome,
      subscale  = sc_raw,
      r         = ct$estimate,
      ci_lower  = ct$conf.int[1],
      ci_upper  = ct$conf.int[2],
      p_value   = ct$p.value,
      n         = nrow(d)
    )
  })
})

# BH correction within each outcome
subscale_cors <- subscale_cors |>
  dplyr::group_by(outcome) |>
  dplyr::mutate(p_BH = p.adjust(p_value, method = "BH"),
                sig_BH = p_BH < .05) |>
  dplyr::ungroup()

cat("\nMAIA subscale correlations with threshold and confidence:\n")
print(as.data.frame(subscale_cors), digits = 3, row.names = FALSE)

readr::write_csv(subscale_cors,
                 file.path(RESULTS_DIR, "s5e_maia_subscale_cors.csv"))
message("Saved: s5e_maia_subscale_cors.csv")

# ── E2b. Alexithymia dissociation ─────────────────────────────
#
# Convergent test of the metacognitive calibration interpretation.
# If MAIA measures calibration (knowing what you know), then
# alexithymia (difficulty identifying and describing feelings)
# should mirror the MAIA dissociation:
#   Alexithymia ~ Awareness   (negative: worse calibration)
#   Alexithymia ~ Threshold   (null: detection sensitivity unaffected)
#
# This replicates the H3 dissociation pattern using a clinically
# relevant individual-difference measure, strengthening the
# metacognitive framing of MAIA.
cat("\n--- E2b: Alexithymia dissociation ---\n")

d_alex <- summary_data |>
  dplyr::filter(!is.na(Alexithymia)) |>
  dplyr::mutate(
    alex_z = scale(Alexithymia)[, 1],
    overall_mean_threshold = rowMeans(
      dplyr::across(dplyr::matches("^thresh_ses1_")), na.rm = TRUE)
  )

cat(sprintf("N with Alexithymia: %d\n", nrow(d_alex)))

# Awareness ~ Alexithymia
alex_aware <- cor.test(d_alex$Awareness_ses1, d_alex$alex_z)
cat(sprintf("Awareness_ses1 ~ Alexithymia_z: %s\n", fmt_cor(alex_aware)))

# Threshold ~ Alexithymia (expect null)
alex_thresh <- cor.test(d_alex$overall_mean_threshold, d_alex$alex_z)
cat(sprintf("Threshold      ~ Alexithymia_z: %s\n", fmt_cor(alex_thresh)))

# Confidence ~ Alexithymia
alex_conf <- cor.test(d_alex$mean_Confidence_ses1, d_alex$alex_z)
cat(sprintf("Confidence     ~ Alexithymia_z: %s\n", fmt_cor(alex_conf)))

# Bayesian null for threshold association
alex_thresh_bf <- tryCatch(
  BayesFactor::correlationBF(d_alex$overall_mean_threshold,
                               d_alex$Alexithymia, iterations = 10000),
  error = function(e) NULL
)
if (!is.null(alex_thresh_bf)) {
  bf10 <- BayesFactor::extractBF(alex_thresh_bf)$bf
  cat(sprintf("Threshold ~ Alexithymia BF01 = %.3f  (%s)\n",
              1/bf10,
              dplyr::case_when(
                1/bf10 >= 10 ~ "strong null evidence",
                1/bf10 >= 3  ~ "moderate null evidence",
                TRUE         ~ "anecdotal")))
}

cat("\nInterpretation:\n")
cat("  If Awareness ~ Alexithymia (negative) but Threshold ~ null:\n")
cat("  → replicates H3 dissociation with clinical measure\n")
cat("  → supports MAIA-as-metacognitive-calibration framing\n")

alex_results <- tibble::tibble(
  outcome     = c("Awareness_ses1", "Threshold", "Confidence_ses1"),
  r           = c(alex_aware$estimate, alex_thresh$estimate,
                  alex_conf$estimate),
  ci_lower    = c(alex_aware$conf.int[1], alex_thresh$conf.int[1],
                  alex_conf$conf.int[1]),
  ci_upper    = c(alex_aware$conf.int[2], alex_thresh$conf.int[2],
                  alex_conf$conf.int[2]),
  p_value     = c(alex_aware$p.value, alex_thresh$p.value,
                  alex_conf$p.value),
  n           = nrow(d_alex)
)
readr::write_csv(alex_results,
                 file.path(RESULTS_DIR, "s5e_alexithymia_dissociation.csv"))
message("Saved: s5e_alexithymia_dissociation.csv")
#
# H3A tests whether MAIA correlates with overall mean Confidence.
# Here we ask whether MAIA correlates with the Awareness index
# (within-person cor(Confidence, Accuracy)) specifically, and
# whether this differs by session or breathing condition.
# ============================================================
cat("\n\n========================================\n")
cat("E3: AWARENESS × MAIA — SESSION/CONDITION BREAKDOWN\n")
cat("========================================\n")

awareness_maia <- summary_data |>
  dplyr::filter(!is.na(MAIA_total_z)) |>
  dplyr::select(id, Group, MAIA_total_z,
                Awareness_ses1, Awareness_ses2)

for (ses_col in c("Awareness_ses1", "Awareness_ses2")) {
  d <- awareness_maia |> dplyr::filter(!is.na(.data[[ses_col]]))
  ct <- cor.test(d[[ses_col]], d$MAIA_total_z)
  cat(sprintf("%s ~ MAIA: %s  N=%d\n",
              ses_col, fmt_cor(ct), nrow(d)))
}

# Cross-session stability of Awareness
d_both <- awareness_maia |>
  dplyr::filter(!is.na(Awareness_ses1), !is.na(Awareness_ses2))
cat(sprintf("\nAwareness ses1 ~ ses2 correlation: %s  N=%d\n",
            fmt_cor(cor.test(d_both$Awareness_ses1, d_both$Awareness_ses2)),
            nrow(d_both)))


# ============================================================
# E4. GROUP DIFFERENCES: BREATH-FIRST VS VISUAL-FIRST
#
# Do Breath-first participants develop higher Awareness by ses2?
# This would suggest that sustained interoceptive attention
# (breathing with the circle) has a training-like effect on
# metacognitive accuracy — consistent with the theoretical claim
# that the task trains interoceptive awareness.
# ============================================================
cat("\n\n========================================\n")
cat("E4: GROUP DIFFERENCES IN AWARENESS\n")
cat("Breath-first vs Visual-first\n")
cat("========================================\n")

group_awareness <- summary_data |>
  dplyr::filter(!is.na(Group)) |>
  dplyr::select(id, Group, Awareness_ses1, Awareness_ses2,
                mean_Confidence_ses1, mean_Confidence_ses2)

# Awareness ses1
cat("\nSession 1 Awareness by Group:\n")
group_awareness |>
  dplyr::filter(!is.na(Awareness_ses1)) |>
  dplyr::group_by(Group) |>
  dplyr::summarise(M = mean(Awareness_ses1),
                   SD = sd(Awareness_ses1), n = dplyr::n(),
                   .groups = "drop") |>
  print()

m_grp1 <- lm(Awareness_ses1 ~ Group, data = group_awareness)
cat(sprintf("Group effect ses1: b=%.3f  p=%.4f\n",
            coef(m_grp1)["Groupvisual"],
            summary(m_grp1)$coefficients["Groupvisual", "Pr(>|t|)"]))

# Awareness ses2 — Visual-first participants doing breath for first time
cat("\nSession 2 Awareness by Group (Visual-first doing breath for 1st time):\n")
group_awareness |>
  dplyr::filter(!is.na(Awareness_ses2)) |>
  dplyr::group_by(Group) |>
  dplyr::summarise(M = mean(Awareness_ses2),
                   SD = sd(Awareness_ses2), n = dplyr::n(),
                   .groups = "drop") |>
  print()

m_grp2 <- lm(Awareness_ses2 ~ Group,
             data = group_awareness |> dplyr::filter(!is.na(Awareness_ses2)))
cat(sprintf("Group effect ses2: b=%.3f  p=%.4f\n",
            coef(m_grp2)["Groupvisual"],
            summary(m_grp2)$coefficients["Groupvisual", "Pr(>|t|)"]))

# Save
group_aware_table <- group_awareness |>
  tidyr::pivot_longer(
    cols = c(Awareness_ses1, Awareness_ses2),
    names_to = "ses", values_to = "Awareness"
  ) |>
  dplyr::filter(!is.na(Awareness)) |>
  dplyr::group_by(Group, ses) |>
  dplyr::summarise(M = mean(Awareness), SD = sd(Awareness),
                   n = dplyr::n(), .groups = "drop")
readr::write_csv(group_aware_table,
                 file.path(RESULTS_DIR, "s5e_group_awareness.csv"))
message("Saved: s5e_group_awareness.csv")

# ── Training contrast: first breath-pacing session ────────────
#
# Tests whether a single session of breath-pacing has a
# metacognitive training effect, independent of the broader
# experimental context.
#
# Key comparison: Awareness at a participant's FIRST breath-pacing
# session regardless of group assignment:
#   Breath-first (Group="breath"): first breath session = ses1
#   Visual-first  (Group="visual"): first breath session = ses2
#
# If these are equivalent → no carry-over from visual condition,
#   and there is no single-session learning effect
# If Breath-first ses1 > Visual-first ses2 → experimental context
#   (having done breath-pacing before, even implicitly) matters
# If Breath-first ses1 ≈ Visual-first ses2, and both > some
#   reference → single breath session produces the same Awareness
#   regardless of what came before

cat("\n--- E4b: Training contrast (first breath-pacing session) ---\n")

first_breath <- group_awareness |>
  dplyr::mutate(
    # First breath session Awareness for each participant
    Awareness_first_breath = dplyr::if_else(
      Group == "breath", Awareness_ses1, Awareness_ses2
    )
  ) |>
  dplyr::filter(!is.na(Awareness_first_breath))

cat(sprintf("N (Breath-first first-breath ses): %d\n",
            sum(first_breath$Group == "breath")))
cat(sprintf("N (Visual-first first-breath ses): %d\n",
            sum(first_breath$Group == "visual")))

first_breath |>
  dplyr::group_by(Group) |>
  dplyr::summarise(
    M_first_breath = mean(Awareness_first_breath),
    SD             = sd(Awareness_first_breath),
    n              = dplyr::n(),
    .groups = "drop"
  ) |>
  print()

m_first_breath <- lm(Awareness_first_breath ~ Group,
                     data = first_breath)
cat(sprintf("\nGroup difference at first breath session: b=%.3f  SE=%.3f  p=%.4f\n",
            coef(m_first_breath)["Groupvisual"],
            summary(m_first_breath)$coefficients["Groupvisual", "Std. Error"],
            summary(m_first_breath)$coefficients["Groupvisual", "Pr(>|t|)"]))
cat("  b ≈ 0 → single breath session produces equivalent Awareness regardless of history\n")
cat("  b > 0 (Breath > Visual at first session) → experimental context matters\n")

# Planned contrast: Breath-first ses1 vs Visual-first ses1
# (equivalent to the pure context effect before any breath-pacing for either)
cat(sprintf("\nContext-only contrast (Visual ses1 = pre-breath baseline):\n"))
context_d <- group_awareness |>
  dplyr::filter(!is.na(Awareness_ses1)) |>
  dplyr::group_by(Group) |>
  dplyr::summarise(M = mean(Awareness_ses1), SD = sd(Awareness_ses1),
                   n = dplyr::n(), .groups = "drop")
print(as.data.frame(context_d), digits = 3, row.names = FALSE)

# Export model summary alongside raw data
e4b_coef <- summary(m_first_breath)$coefficients
e4b_results <- tibble::tibble(
  analysis        = "E4b: Training contrast (first breath session)",
  b_Groupvisual   = e4b_coef["Groupvisual", "Estimate"],
  se_Groupvisual  = e4b_coef["Groupvisual", "Std. Error"],
  t_Groupvisual   = e4b_coef["Groupvisual", "t value"],
  p_Groupvisual   = e4b_coef["Groupvisual", "Pr(>|t|)"],
  n_breath        = sum(first_breath$Group == "breath"),
  n_visual        = sum(first_breath$Group == "visual"),
  M_breath        = mean(first_breath$Awareness_first_breath[
                      first_breath$Group == "breath"], na.rm = TRUE),
  M_visual        = mean(first_breath$Awareness_first_breath[
                      first_breath$Group == "visual"], na.rm = TRUE)
)
readr::write_csv(e4b_results,
  file.path(RESULTS_DIR, "s5e_training_contrast.csv")
)
message("Saved: s5e_training_contrast.csv")
#
# Does Awareness improve from ses1 to ses2?
# If so, is the improvement larger for Visual-first participants
# whose ses2 is their first breath-pacing session?
# ============================================================
cat("\n\n========================================\n")
cat("E5: AWARENESS CHANGE SESSION 1 → SESSION 2\n")
cat("========================================\n")

d_long_aware <- group_awareness |>
  dplyr::filter(!is.na(Awareness_ses1), !is.na(Awareness_ses2)) |>
  dplyr::mutate(Awareness_change = Awareness_ses2 - Awareness_ses1)

cat(sprintf(
  "N with both sessions: %d  |  Mean change: %.3f (SD=%.3f)\n",
  nrow(d_long_aware),
  mean(d_long_aware$Awareness_change),
  sd(d_long_aware$Awareness_change)
))

# Overall change
tt_change <- t.test(d_long_aware$Awareness_ses2,
                    d_long_aware$Awareness_ses1,
                    paired = TRUE)
cat(sprintf("Paired t-test ses1 vs ses2: t(%.0f)=%.3f  p=%.4f\n",
            tt_change$parameter, tt_change$statistic, tt_change$p.value))

# By group: is the change larger for Visual-first?
m_change <- lm(Awareness_change ~ Group, data = d_long_aware)
cat(sprintf("Group × change interaction: b=%.3f  p=%.4f\n",
            coef(m_change)["Groupvisual"],
            summary(m_change)$coefficients["Groupvisual", "Pr(>|t|)"]))

d_long_aware |>
  dplyr::group_by(Group) |>
  dplyr::summarise(
    M_change = mean(Awareness_change),
    SD_change = sd(Awareness_change),
    n = dplyr::n(), .groups = "drop"
  ) |>
  print()

# Export model summary alongside raw data
e5_change_coef <- summary(m_change)$coefficients
e5_results <- tibble::tibble(
  analysis           = "E5: Awareness change ses1 to ses2",
  n_paired           = nrow(d_long_aware),
  mean_change        = mean(d_long_aware$Awareness_change),
  sd_change          = sd(d_long_aware$Awareness_change),
  t_paired           = unname(tt_change$statistic),
  df_paired          = unname(tt_change$parameter),
  p_paired           = tt_change$p.value,
  b_Groupvisual      = e5_change_coef["Groupvisual", "Estimate"],
  se_Groupvisual     = e5_change_coef["Groupvisual", "Std. Error"],
  p_Groupvisual      = e5_change_coef["Groupvisual", "Pr(>|t|)"],
  M_change_breath    = mean(d_long_aware$Awareness_change[
                         d_long_aware$Group == "breath"], na.rm = TRUE),
  M_change_visual    = mean(d_long_aware$Awareness_change[
                         d_long_aware$Group == "visual"], na.rm = TRUE)
)
readr::write_csv(e5_results,
  file.path(RESULTS_DIR, "s5e_awareness_change.csv")
)
message("Saved: s5e_awareness_change.csv")


# ============================================================
# FIGURES
# ============================================================

# F1: H3B breakdown — r across conditions
f_e1 <- dplyr::bind_rows(
  tibble::tibble(condition = "Overall", r = e1_overall$estimate,
                 ci_lower  = e1_overall$conf.int[1],
                 ci_upper  = e1_overall$conf.int[2],
                 level     = "Primary"),
  level2_results |>
    dplyr::mutate(condition = paste(factor, level), level = "Simple effects") |>
    dplyr::select(condition, r, ci_lower, ci_upper, level),
  level3_results |>
    dplyr::mutate(
      condition = stringr::str_remove(condition, "thresh_ses1_"),
      level     = "Cells"
    ) |>
    dplyr::select(condition, r, ci_lower, ci_upper, level)
) |>
  dplyr::mutate(
    level     = factor(level, levels = c("Primary", "Simple effects", "Cells")),
    condition = factor(condition, levels = rev(unique(condition)))
  ) |>
  ggplot(aes(x = r, y = condition,
             xmin = ci_lower, xmax = ci_upper,
             colour = level)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(height = 0.2, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "D", begin = 0.1, end = 0.8) +
  labs(
    title    = "H3B: Threshold ~ MAIA across conditions",
    subtitle = "Exploratory breakdown — all conditions expected null",
    x = "r (threshold ~ MAIA_z)",
    y = NULL, colour = "Level"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(s5e_fig_dir, "s5e_h3b_breakdown.pdf"),
       f_e1, width = 8, height = 6, device = "pdf")
message("Saved: s5e_h3b_breakdown.pdf")

# F2: Awareness change by group
f_e5 <- d_long_aware |>
  tidyr::pivot_longer(c(Awareness_ses1, Awareness_ses2),
                      names_to = "ses", values_to = "Awareness") |>
  dplyr::mutate(
    ses = factor(ses, levels = c("Awareness_ses1", "Awareness_ses2"),
                 labels = c("Session 1", "Session 2")),
    Group = factor(Group, levels = c("breath", "visual"),
                   labels = c("Breath-first", "Visual-first"))
  ) |>
  ggplot(aes(x = ses, y = Awareness, group = id,
             colour = Group)) +
  geom_line(alpha = 0.25, linewidth = 0.4) +
  stat_summary(aes(group = Group),
               fun = mean, geom = "line",
               linewidth = 1.5, alpha = 0.9) +
  stat_summary(aes(group = Group),
               fun.data = mean_se, geom = "pointrange",
               linewidth = 1.0, size = 0.6) +
  scale_colour_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  labs(
    title    = "Awareness index: session 1 → session 2",
    subtitle = "Thin lines = individual participants; thick = group mean ± SE",
    x = NULL, y = "Awareness (r)",
    colour = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(s5e_fig_dir, "s5e_awareness_change.pdf"),
       f_e5, width = 7, height = 5, device = "pdf")
message("Saved: s5e_awareness_change.pdf")

message("\nanalysis_study5_exploratory.R complete.")
message("Results written to: ", RESULTS_DIR)
