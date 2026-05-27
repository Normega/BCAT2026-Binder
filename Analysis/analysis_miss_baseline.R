# ============================================================
# analysis_miss_baseline.R
# Bayesian multilevel null test: do missed-change trials produce
# arousal above the no-change baseline?
#
# Runs across all four studies with no-change baselines:
#   Studies 1A, 2, 4, 5
# (Study 3 excluded: no true no-change condition)
#
# Approach: for each study, fit two brms models on
# {no-change trials + missed change trials}:
#   M0: Arousal ~ 1           + (1 | id)   [intercept only]
#   M1: Arousal ~ IsMiss      + (1 | id)   [miss vs. baseline]
# Bayes factor via bridge sampling (BF10 = M1 / M0).
# BF01 > 3 = moderate evidence for null (no elevation on misses).
#
# Person-level ttestBF also reported for comparison with
# previously reported values (Studies 4/5: BF01 = 4.61, 5.41).
#
# Assumes in environment (from MainAnalysis.R):
#   s1l, s2l, s4l, s5_long_breath
#   RESULTS_DIR
#
# Output:
#   miss_baseline_bf.csv — BF01 null tests per study (multilevel + person-level)
# ============================================================


message("\n========================================")
message("MISS vs. BASELINE: Multilevel Bayesian Null Test")
message("========================================")


# ── Helper: prep miss + no-change data for one study ─────────
.prep_miss_baseline <- function(long_data, study_label) {

  dat <- long_data |>
    dplyr::filter(
      !is.na(Arousal),
      !is.na(Accuracy),
      Direction %in% c("NoChange", "Faster", "Slower")
    ) |>
    dplyr::mutate(
      IsMiss = dplyr::case_when(
        Direction == "NoChange"              ~ 0L,  # no-change baseline
        Direction != "NoChange" & Accuracy == 0 ~ 1L,  # missed change
        TRUE ~ NA_integer_
      )
    ) |>
    dplyr::filter(!is.na(IsMiss))

  n_nochange <- sum(dat$IsMiss == 0)
  n_miss     <- sum(dat$IsMiss == 1)
  n_subj     <- dplyr::n_distinct(dat$id)

  message(sprintf("\n[%s] No-change trials: %d | Miss trials: %d | N participants: %d",
                  study_label, n_nochange, n_miss, n_subj))

  dat
}


# ── Helper: run multilevel BF test for one study ──────────────
.run_miss_bf <- function(dat, study_label,
                          chains = 4, iter = 4000, warmup = 1000,
                          seed = 42) {

  # Scale arousal within-study for prior calibration
  dat <- dat |>
    dplyr::mutate(Arousal_z = as.numeric(scale(Arousal)),
                  IsMiss    = as.numeric(IsMiss))

  # Priors: JZS-style Cauchy on fixed effect (r = 0.707),
  # half-Cauchy on random effect SD and residual SD
  priors <- c(
    brms::prior(cauchy(0, 0.707), class = b),
    brms::prior(cauchy(0, 1),     class = sd),
    brms::prior(cauchy(0, 1),     class = sigma)
  )

  message(sprintf("  Fitting M0 (intercept only) for %s...", study_label))
  m0 <- brms::brm(
    Arousal_z ~ 1 + (1 | id),
    data      = dat,
    prior     = priors[priors$class != "b", ],
    chains    = chains,
    iter      = iter,
    warmup    = warmup,
    seed      = seed,
    save_pars = brms::save_pars(all = TRUE),
    silent    = 2
  )

  message(sprintf("  Fitting M1 (IsMiss) for %s...", study_label))
  m1 <- brms::brm(
    Arousal_z ~ IsMiss + (1 | id),
    data      = dat,
    prior     = priors,
    chains    = chains,
    iter      = iter,
    warmup    = warmup,
    seed      = seed,
    save_pars = brms::save_pars(all = TRUE),
    silent    = 2
  )

  message(sprintf("  Bridge sampling for %s...", study_label))
  bf_obj <- brms::bayes_factor(m1, m0)
  bf10   <- bf_obj$bf
  bf01   <- 1 / bf10

  # Frequentist LME for reference
  m_freq <- lmerTest::lmer(Arousal ~ IsMiss + (1 | id),
                            data = dat, REML = FALSE)
  cf     <- summary(m_freq)$coefficients
  b_miss <- cf["IsMiss", "Estimate"]
  se_miss <- cf["IsMiss", "Std. Error"]
  t_miss  <- cf["IsMiss", "t value"]
  p_miss  <- cf["IsMiss", "Pr(>|t|)"]

  # Person-level ttestBF for comparison with paper's existing values
  person_diffs <- dat |>
    dplyr::group_by(id) |>
    dplyr::summarise(
      mean_miss     = mean(Arousal[IsMiss == 1], na.rm = TRUE),
      mean_nochange = mean(Arousal[IsMiss == 0], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::filter(!is.na(mean_miss), !is.na(mean_nochange)) |>
    dplyr::mutate(diff = mean_miss - mean_nochange)

  bf_person     <- BayesFactor::ttestBF(person_diffs$diff, mu = 0)
  bf10_person   <- exp(bf_person@bayesFactor$bf)
  bf01_person   <- 1 / bf10_person

  tibble::tibble(
    study        = study_label,
    n_nochange   = sum(dat$IsMiss == 0),
    n_miss       = sum(dat$IsMiss == 1),
    n_subj       = dplyr::n_distinct(dat$id),
    b_miss       = b_miss,
    se_miss      = se_miss,
    t_miss       = t_miss,
    p_miss       = p_miss,
    BF10_multilevel = bf10,
    BF01_multilevel = bf01,
    BF10_person  = bf10_person,
    BF01_person  = bf01_person,
    n_person_bf  = nrow(person_diffs)
  )
}


# ── Prep datasets ─────────────────────────────────────────────
dat_s1a <- .prep_miss_baseline(
  dplyr::filter(s1l, Group == "TaskA"), "Study1A")
dat_s2  <- .prep_miss_baseline(s2l,            "Study2")
dat_s4  <- .prep_miss_baseline(s4l,            "Study4")
dat_s5  <- .prep_miss_baseline(s5_long_breath, "Study5")


# ── Run BF tests ──────────────────────────────────────────────
# Note: brms sampling takes several minutes per study.
# Set chains = 2, iter = 2000 for a quick check run.

results <- purrr::map_dfr(
  list(
    list(dat = dat_s1a, label = "Study1A"),
    list(dat = dat_s2,  label = "Study2"),
    list(dat = dat_s4,  label = "Study4"),
    list(dat = dat_s5,  label = "Study5")
  ),
  function(x) {
    tryCatch(
      .run_miss_bf(x$dat, x$label),
      error = function(e) {
        message(sprintf("  ERROR in %s: %s", x$label, e$message))
        tibble::tibble(study = x$label)
      }
    )
  }
)


# ── Print and save ────────────────────────────────────────────
message("\n--- Results: Miss vs. No-Change Arousal ---")
print(
  results |>
    dplyr::mutate(
      across(c(b_miss, se_miss, t_miss, p_miss), \(x) round(x, 3)),
      across(c(BF10_multilevel, BF01_multilevel,
               BF10_person,    BF01_person),    \(x) round(x, 2))
    )
)

readr::write_csv(results,
                 file.path(RESULTS_DIR, "miss_baseline_bf.csv"))
message("Saved: miss_baseline_bf.csv")


# ── Interpretation guide ──────────────────────────────────────
message("\n--- Interpretation ---")
message("BF01 > 3  = moderate evidence for null (no arousal elevation on misses)")
message("BF01 1-3  = inconclusive")
message("BF01 < 1  = evidence for elevation on misses")
message("")
message("Primary: BF01_multilevel (consistent with all other arousal models)")
message("Reference: BF01_person (matches previously reported Study 4/5 values)")
