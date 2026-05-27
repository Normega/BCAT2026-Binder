# ============================================================
# analysis_hbd.R
# HBD cardiac interoception analyses (Study 5, exploratory).
# Sources after analysis_maia.R; uses s5s, s5_long_breath.
# Outputs: hbd_missingness.csv, hbd3_crossmodal.csv,
#          hbd5_sensitivity.csv, hbd_summary.csv
# ============================================================

source(file.path(ANALYSIS_DIR, "theme_bcat.R"))

hbd_fig_dir <- file.path(FIG_DIR, "Study5_HBD")
dir.create(hbd_fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load ──────────────────────────────────────────────────────
# s5s and s5_long_breath from MainAnalysis.R
summary_data <- s5s
long_data    <- s5_long_breath

# HBD data loaded by MainAnalysis.R (via utils.R load_all_data)
hbd_sum       <- s5_hbd
hbd_intervals <- s5_hbd_intervals

# Sensitivity Schandry variants
schandry_sens <- hbd_intervals |>
  dplyr::group_by(id) |>
  dplyr::summarise(
    schandry_25s    = mean(schandry_score[interval_num == 1 &
                                           hr_plausible != FALSE], na.rm=TRUE),
    schandry_25_35s = mean(schandry_score[interval_num %in% c(1,2) &
                                           hr_plausible != FALSE], na.rm=TRUE),
    .groups = "drop")

bcat_sum <- summary_data |>
  dplyr::mutate(
    overall_mean_threshold = rowMeans(
      dplyr::across(dplyr::matches("^thresh_ses1_")), na.rm = TRUE),
    mean_confidence = mean_Confidence_ses1
  ) |>
  dplyr::select(id, Group, overall_mean_threshold, mean_confidence,
                Awareness_ses1, MAIA_total, MAIA_total_z,
                MSES_selfdoubt, Alexithymia)

merged <- bcat_sum |>
  dplyr::left_join(
    hbd_sum |> dplyr::select(id, mean_schandry, mean_abs_error,
                               mean_hr_bpm, onset_source,
                               n_intervals_analyzed, hbd_quality,
                               hbd_implausible),
    by = "id") |>
  dplyr::left_join(schandry_sens, by = "id") |>
  dplyr::filter(!hbd_implausible)

message(sprintf("Merged: N=%d | complete HBD=%d (%.0f%%)",
  nrow(merged),
  sum(merged$hbd_quality=="complete", na.rm=TRUE),
  100*mean(merged$hbd_quality=="complete", na.rm=TRUE)))

# ── HBD-0. Missingness ───────────────────────────────────────
cat("\n=== HBD-0: Missingness analysis ===\n")
miss_vars <- c("overall_mean_threshold","mean_confidence",
               "MAIA_total_z","Awareness_ses1")
miss_results <- purrr::map_dfr(miss_vars, function(v) {
  if (!v %in% names(merged)) return(NULL)
  present <- merged |> dplyr::filter(hbd_quality=="complete",
                                      !is.na(.data[[v]])) |> dplyr::pull(v)
  absent  <- merged |> dplyr::filter(hbd_quality!="complete",
                                      !is.na(.data[[v]])) |> dplyr::pull(v)
  if (length(present)<3 || length(absent)<3) return(NULL)
  tt <- t.test(present, absent)
  tibble::tibble(variable=v, M_present=mean(present), M_absent=mean(absent),
    t=tt$statistic, df=tt$parameter, p=tt$p.value,
    d=(mean(present)-mean(absent))/sqrt((var(present)+var(absent))/2))
})
print(as.data.frame(miss_results), digits=3, row.names=FALSE)
readr::write_csv(miss_results, file.path(RESULTS_DIR, "hbd_missingness.csv"))

# ── HBD-1. Convergent validity ───────────────────────────────
cat("\n=== HBD-1: Schandry ~ BCAT threshold ===\n")
d1 <- merged |> dplyr::filter(!is.na(mean_schandry),
                                !is.na(overall_mean_threshold))
cat(sprintf("N=%d\n", nrow(d1)))
hbd1_zero    <- cor.test(d1$mean_schandry, d1$overall_mean_threshold)
cat(sprintf("Zero-order: %s\n", fmt_cor(hbd1_zero)))
d1p <- d1 |> dplyr::filter(!is.na(mean_hr_bpm))
hbd1_partial <- cor.test(
  residuals(lm(mean_schandry~mean_hr_bpm, data=d1p)),
  residuals(lm(overall_mean_threshold~mean_hr_bpm, data=d1p)))
cat(sprintf("Partial (|HR): %s  N=%d\n", fmt_cor(hbd1_partial), nrow(d1p)))
hbd1_reg <- lm(mean_schandry ~ overall_mean_threshold + mean_hr_bpm,
               data=d1p)
print(summary(hbd1_reg))

# ── HBD-2. MAIA dissociation at cardiac level ─────────────────
cat("\n=== HBD-2: Schandry ~ MAIA (expect null) ===\n")
d2 <- merged |> dplyr::filter(!is.na(mean_schandry), !is.na(MAIA_total_z))
cat(sprintf("N=%d\n", nrow(d2)))
hbd2_zero <- cor.test(d2$mean_schandry, d2$MAIA_total_z)
cat(sprintf("Zero-order: %s\n", fmt_cor(hbd2_zero)))
hbd2_bf <- tryCatch(
  BayesFactor::correlationBF(d2$mean_schandry, d2$MAIA_total,
                               iterations=10000),
  error=function(e) NULL)
if (!is.null(hbd2_bf)) {
  bf10 <- BayesFactor::extractBF(hbd2_bf)$bf
  cat(sprintf("BF10=%.3f  BF01=%.3f\n", bf10, 1/bf10))
}
# Interval-level model
d2_int <- hbd_intervals |>
  dplyr::left_join(d2|>dplyr::select(id,MAIA_total_z), by="id") |>
  dplyr::filter(!is.na(schandry_score),!is.na(MAIA_total_z),
                hr_plausible!=FALSE)
m_hbd2_int <- lmerTest::lmer(
  schandry_score ~ MAIA_total_z + factor(interval_num) + (1|id),
  data=d2_int, REML=FALSE)
cat("\nInterval-level model:\n")
print(broom.mixed::tidy(m_hbd2_int, effects="fixed", conf.int=TRUE), n=Inf)

# ── HBD-3. Cross-modal confidence (primary novel contribution) ──
cat("\n=== HBD-3: Cross-modal confidence (primary novel) ===\n")
d3 <- merged |> dplyr::filter(!is.na(mean_schandry),
                                !is.na(mean_confidence),
                                !is.na(overall_mean_threshold))
cat(sprintf("N=%d\n", nrow(d3)))
hbd3_zero    <- cor.test(d3$mean_schandry, d3$mean_confidence)
cat(sprintf("Zero-order confidence~Schandry: %s\n", fmt_cor(hbd3_zero)))
hbd3_partial <- cor.test(
  residuals(lm(mean_confidence~overall_mean_threshold, data=d3)),
  residuals(lm(mean_schandry~overall_mean_threshold,   data=d3)))
cat(sprintf("Partial (|threshold): %s\n", fmt_cor(hbd3_partial)))
hbd3_reg <- lm(mean_schandry ~ mean_confidence + overall_mean_threshold,
               data=d3)
cat("\nRegression (Schandry ~ confidence + threshold):\n")
print(summary(hbd3_reg))
if ("MAIA_total_z" %in% names(d3)) {
  d3m <- d3 |> dplyr::filter(!is.na(MAIA_total_z))
  cat("\nRobustness (+ MAIA):\n")
  print(summary(lm(mean_schandry ~ mean_confidence + overall_mean_threshold +
                     MAIA_total_z, data=d3m)))
}
hbd3_table <- tibble::tibble(
  analysis = c("Zero-order","Partial (|threshold)","Regression β_confidence"),
  estimate = c(hbd3_zero$estimate, hbd3_partial$estimate,
               coef(hbd3_reg)["mean_confidence"]),
  ci_lower = c(hbd3_zero$conf.int[1], hbd3_partial$conf.int[1],
               confint(hbd3_reg)["mean_confidence",1]),
  ci_upper = c(hbd3_zero$conf.int[2], hbd3_partial$conf.int[2],
               confint(hbd3_reg)["mean_confidence",2]),
  p_value  = c(hbd3_zero$p.value, hbd3_partial$p.value,
               summary(hbd3_reg)$coefficients["mean_confidence","Pr(>|t|)"]))
readr::write_csv(hbd3_table, file.path(RESULTS_DIR, "hbd3_crossmodal.csv"))

# ── HBD-4. Arousal moderation (exploratory) ──────────────────
cat("\n=== HBD-4: Arousal ~ Change × Schandry (exploratory) ===\n")
long_h4 <- long_data |>
  dplyr::filter(Condition=="breath", !is.na(Arousal), !is.na(Change)) |>
  dplyr::left_join(merged|>dplyr::select(id,mean_schandry), by="id") |>
  dplyr::filter(!is.na(mean_schandry)) |>
  dplyr::mutate(Change2=Change^2,
                schandry_z=scale(mean_schandry)[,1])
cat(sprintf("N=%d trials, %d participants\n",
            nrow(long_h4), dplyr::n_distinct(long_h4$id)))
m_h4_base <- lmerTest::lmer(Arousal~Change+Change2+(1|id), data=long_h4, REML=FALSE)
m_h4_mod  <- lmerTest::lmer(Arousal~Change*schandry_z+Change2+(1|id),
                              data=long_h4, REML=FALSE)
cat("\nLRT:\n"); print(anova(m_h4_base, m_h4_mod))
cat("\nFull model:\n")
print(broom.mixed::tidy(m_h4_mod, effects="fixed", conf.int=TRUE), n=Inf)

# ── HBD-5. Sensitivity ───────────────────────────────────────
cat("\n=== HBD-5: Sensitivity (interval subsets) ===\n")
sens_table <- purrr::map_dfr(
  c("mean_schandry","schandry_25s","schandry_25_35s"),
  function(sv) {
    ds <- merged |> dplyr::filter(!is.na(.data[[sv]]),
                                   !is.na(overall_mean_threshold),
                                   !is.na(mean_confidence))
    rt <- cor.test(ds[[sv]], ds$overall_mean_threshold)
    rc <- cor.test(ds[[sv]], ds$mean_confidence)
    rm <- if("MAIA_total_z" %in% names(ds))
            cor.test(ds[[sv]], ds$MAIA_total_z) else NULL
    tibble::tibble(
      schandry_var = sv, n = nrow(ds),
      r_threshold  = rt$estimate, p_threshold = rt$p.value,
      r_confidence = rc$estimate, p_confidence = rc$p.value,
      r_maia       = if(!is.null(rm)) rm$estimate else NA,
      p_maia       = if(!is.null(rm)) rm$p.value  else NA)
  })
print(as.data.frame(sens_table), digits=3, row.names=FALSE)
readr::write_csv(sens_table, file.path(RESULTS_DIR, "hbd5_sensitivity.csv"))

# ── Summary CSV ──────────────────────────────────────────────
hbd_summary <- tibble::tibble(
  analysis   = c("HBD-1 zero-order","HBD-1 partial|HR",
                 "HBD-2 MAIA zero-order",
                 "HBD-3 confidence zero-order",
                 "HBD-3 confidence partial|threshold"),
  r          = c(hbd1_zero$estimate, hbd1_partial$estimate,
                 hbd2_zero$estimate, hbd3_zero$estimate,
                 hbd3_partial$estimate),
  ci_lower   = c(hbd1_zero$conf.int[1], hbd1_partial$conf.int[1],
                 hbd2_zero$conf.int[1], hbd3_zero$conf.int[1],
                 hbd3_partial$conf.int[1]),
  ci_upper   = c(hbd1_zero$conf.int[2], hbd1_partial$conf.int[2],
                 hbd2_zero$conf.int[2], hbd3_zero$conf.int[2],
                 hbd3_partial$conf.int[2]),
  p_value    = c(hbd1_zero$p.value, hbd1_partial$p.value,
                 hbd2_zero$p.value, hbd3_zero$p.value,
                 hbd3_partial$p.value),
  n          = c(nrow(d1), nrow(d1p), nrow(d2), nrow(d3), nrow(d3)))
readr::write_csv(hbd_summary, file.path(RESULTS_DIR, "hbd_summary.csv"))

# ── Figure ────────────────────────────────────────────────────
f_hbd <- hbd_summary |>
  dplyr::mutate(analysis=factor(analysis,levels=rev(analysis)),
                sig=p_value<.05) |>
  ggplot(aes(x=r, y=analysis, xmin=ci_lower, xmax=ci_upper, colour=sig)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey50") +
  geom_errorbarh(height=0.2, linewidth=0.8) +
  geom_point(size=3) +
  scale_colour_manual(values=c(`FALSE`="grey60",`TRUE`="#2c7be5"),
                      labels=c(`FALSE`="n.s.",`TRUE`="p<.05"), name=NULL) +
  labs(title="HBD: Cardiac Schandry correlations",
       subtitle="r ± 95% CI  |  Study 5 in-person",
       x="Correlation (r)", y=NULL) +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold"), legend.position="bottom")
ggsave(file.path(hbd_fig_dir, "hbd_summary_forest.pdf"),
       f_hbd, width=9, height=6, device="pdf")
message("analysis_hbd.R complete.")
