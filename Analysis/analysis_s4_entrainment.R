# ============================================================
# analysis_s4_entrainment.R
# Supplementary S4: Belt Entrainment and Interoceptive
# Contribution to Detection.
#
# Sources after analysis_belt.R.
# Uses qcFull and s5s already loaded by MainAnalysis.R.
# Output: s4_entrainment_results.csv
# ============================================================

message("\nS4: Belt Entrainment -> Detection Advantage")

# Data (loaded by MainAnalysis.R) --------------------------------------------
qc      <- qcFull |> dplyr::mutate(id = as.integer(id))
summary <- s5s

# Define breath condition labels (run 1 staircase + run 2 named conditions)
BREATH_CONDITIONS <- c("breath", "lowSalienceAcc", "lowSalienceDec",
                       "highSalienceAcc", "highSalienceDec")

# Compute Per-Person Belt Entrainment Metrics ---------------------------------
# Restrict to: breath conditions, non-unusable belt quality, valid trials only
belt_person <- qc |>
  dplyr::filter(
    condition    %in% BREATH_CONDITIONS,
    belt_quality != "unusable",
    trial_available == TRUE
  ) |>
  dplyr::group_by(id) |>
  dplyr::summarise(
    mean_belt_corr    = mean(correlation,       na.rm = TRUE),
    mean_dir_correct  = mean(direction_correct, na.rm = TRUE),
    mean_mae          = mean(mae,               na.rm = TRUE),
    n_belt_trials     = dplyr::n(),
    .groups = "drop"
  )

cat(sprintf("Participants with belt data: %d\n", nrow(belt_person)))

# =============================================================================
# Analysis 1: Visual-First Group
# Does entrainment predict Breath-over-Visual detection advantage?
# =============================================================================

vf <- summary |>
  dplyr::filter(Condition_ses1 == "visual") |>
  dplyr::mutate(breath_advantage = dprime_ses2 - dprime_ses1)

cat(sprintf("\nVisual-first group: N = %d total, %d with ses2 d'\n",
            nrow(vf),
            sum(!is.na(vf$dprime_ses2))))

# Check dropout: do non-completers differ from completers?
cat("\n--- Dropout check (visual-first) ---\n")
dropout_check <- vf |>
  dplyr::group_by(completed_ses2) |>
  dplyr::summarise(
    n          = dplyr::n(),
    dprime_M   = mean(dprime_ses1, na.rm = TRUE),
    dprime_SD  = sd(dprime_ses1,   na.rm = TRUE),
    MAIA_M     = mean(MAIA_total,  na.rm = TRUE),
    .groups = "drop"
  )
print(dropout_check)

t_drop_dprime <- t.test(dprime_ses1 ~ completed_ses2, data = vf)
t_drop_maia   <- t.test(MAIA_total  ~ completed_ses2, data = vf)
cat(sprintf("dprime_ses1: t(%0.1f) = %.3f, p = %.3f\n",
            t_drop_dprime$parameter, t_drop_dprime$statistic, t_drop_dprime$p.value))
cat(sprintf("MAIA_total:  t(%0.1f) = %.3f, p = %.3f\n",
            t_drop_maia$parameter, t_drop_maia$statistic, t_drop_maia$p.value))

# Merge with belt data
vf_belt <- vf |>
  dplyr::inner_join(belt_person, by = "id") |>
  dplyr::filter(!is.na(breath_advantage))

cat(sprintf("\nVisual-first analysis N (belt + ses2 d'): %d\n", nrow(vf_belt)))

# Correlations: entrainment metrics vs breath advantage
cat("\n--- Visual-first: Entrainment vs Breath-over-Visual advantage ---\n")
metrics_vf <- list(
  mean_belt_corr   = "Belt-pacer correlation",
  mean_dir_correct = "Direction correct",
  mean_mae         = "Mean absolute error"
)

vf_results <- purrr::map_dfr(names(metrics_vf), function(m) {
  dat <- vf_belt |> dplyr::filter(!is.na(.data[[m]]))
  ct  <- cor.test(dat[[m]], dat$breath_advantage)
  tibble::tibble(
    metric  = metrics_vf[[m]],
    r       = ct$estimate,
    t_stat  = ct$statistic,
    df      = ct$parameter,
    p       = ct$p.value,
    n       = nrow(dat),
    ci_lo   = ct$conf.int[1],
    ci_hi   = ct$conf.int[2]
  )
})
print(vf_results)

# =============================================================================
# Analysis 2: Breath-First Group
# Does entrainment predict raw d' within the breath-only session?
# (Reverse causation check: effect should not be inflated vs Analysis 1)
# =============================================================================

bf <- summary |>
  dplyr::filter(Condition_ses1 == "breath")

cat(sprintf("\nBreath-first group: N = %d\n", nrow(bf)))

bf_belt <- bf |>
  dplyr::inner_join(belt_person, by = "id") |>
  dplyr::filter(!is.na(dprime_ses1))

cat(sprintf("Breath-first analysis N (belt + ses1 d'): %d\n", nrow(bf_belt)))

cat("\n--- Breath-first: Entrainment vs raw d' (ses1) ---\n")
bf_results <- purrr::map_dfr(names(metrics_vf), function(m) {
  dat <- bf_belt |> dplyr::filter(!is.na(.data[[m]]))
  ct  <- cor.test(dat[[m]], dat$dprime_ses1)
  tibble::tibble(
    metric  = metrics_vf[[m]],
    r       = ct$estimate,
    t_stat  = ct$statistic,
    df      = ct$parameter,
    p       = ct$p.value,
    n       = nrow(dat),
    ci_lo   = ct$conf.int[1],
    ci_hi   = ct$conf.int[2]
  )
})
print(bf_results)

# =============================================================================
# Summary comparison: Visual-first vs Breath-first (belt-pacer correlation only)
# =============================================================================

cat("\n--- Summary: Belt-pacer correlation effect size comparison ---\n")
cat(sprintf(
  "Visual-first (Breath > Visual advantage): r = %.3f, p = %.3f, N = %d\n",
  vf_results$r[vf_results$metric == "Belt-pacer correlation"],
  vf_results$p[vf_results$metric == "Belt-pacer correlation"],
  vf_results$n[vf_results$metric == "Belt-pacer correlation"]
))
cat(sprintf(
  "Breath-first (raw d'):                    r = %.3f, p = %.3f, N = %d\n",
  bf_results$r[bf_results$metric == "Belt-pacer correlation"],
  bf_results$p[bf_results$metric == "Belt-pacer correlation"],
  bf_results$n[bf_results$metric == "Belt-pacer correlation"]
))
cat("(Similar magnitudes suggest reverse causation is not inflating the visual-first result.)\n")

# Save results ----------------------------------------------------------------
s4_results <- dplyr::bind_rows(
  vf_results |> dplyr::mutate(group = "visual-first"),
  bf_results |> dplyr::mutate(group = "breath-first")
) |>
  dplyr::select(group, metric, r, ci_lo, ci_hi, t_stat, df, p, n)

readr::write_csv(s4_results,
                 file.path(RESULTS_DIR, "s4_entrainment_results.csv"))
message("Saved: s4_entrainment_results.csv")
message("analysis_s4_entrainment.R complete.")
