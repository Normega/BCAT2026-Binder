# ============================================================
# analysis_study3_attraction_mediation.R
# Study 3 misattribution chain: Change -> Arousal -> Attraction.
# Sources after analysis_arousal.R; uses s3l, s3s from MainAnalysis.R.
#
# M0: primary mediation
# M1: + Accuracy moderator (awareness-gating test)
# M2: + fatigue covariates (sensitivity)
# M3: + Salience (sensitivity)
# Frequentist: mediation package. Bayesian: brms.
#
# Outputs:
#   study3_mediation_primary.csv        — primary mediation (M0)
#   study3_mediation_m1a_subset.csv     — M1a accuracy moderation (hit subset)
#   study3_mediation_m1b_analytic.csv   — M1b analytic moderation
#   study3_mediation_bayesian_m0.csv    — Bayesian M0 indirect effects
#   study3_mediation_bayesian_m1.csv    — Bayesian M1 moderated mediation
#   study3_maia_moderation_path_a.csv   — MAIA moderation of path a
#   study3_maia_protection_acme.csv     — MAIA protection: ACME by group
#   study3_maia_bayesian_moderation.csv — Bayesian MAIA moderation
#   Figures/Study3/mediation_forest.pdf — mediation forest plot
# ============================================================

source(file.path(ANALYSIS_DIR, "theme_bcat.R"))

fig_dir   <- file.path(FIG_DIR, "Study3")
cache_dir <- file.path(MODEL_DIR, "Study3")
dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# brms and mediation are loaded by MainAnalysis.R

# ============================================================
# 1. Load and Prepare Data
# ============================================================
# s3l loaded by MainAnalysis.R
s3 <- s3l |>
  dplyr::filter(exclusion_reason == "included") |>
  dplyr::mutate(
    # Quadratic change term (centred; Change already centred at 0)
    Change2      = Change^2,
    # Binary accuracy as numeric for mediation models
    Accuracy     = as.numeric(Accuracy),
    # Sequential within-participant trial counter for fatigue models
    # BlockNum indexes trial within run; Run distinguishes the two halves
    Run_num      = dplyr::if_else(Run == "run1", 0L, 1L),  # 0/1 for model
    # Person-mean centre Change for random-effect interpretability
    # (not needed for M0 but useful for brms random slopes if added later)
    id           = as.integer(id)
  )

message("Participants: ", n_distinct(s3$id))
message("Trials:       ", nrow(s3))
message("Change range: ", round(min(s3$Change), 2), " to ", round(max(s3$Change), 2))

# Contrast values for ACME in mediate()
change_sd  <- sd(s3$Change)
treat_val  <- mean(s3$Change) + change_sd   # +1 SD  (slower breathing)
ctrl_val   <- mean(s3$Change) - change_sd   # −1 SD  (faster breathing)
message(sprintf("ACME contrasts: control=%.3f, treat=%.3f (2-SD window)", ctrl_val, treat_val))

# ============================================================
# 2. Helper: extract and print mediation summary cleanly
# ============================================================
print_med <- function(med_obj, label) {
  s <- summary(med_obj)
  cat("\n──", label, "──\n")
  cat(sprintf("  ACME (indirect):  est=%.4f  95CI [%.4f, %.4f]  p=%.3f\n",
              s$d0, s$d0.ci[1], s$d0.ci[2], s$d0.p))
  cat(sprintf("  ADE  (direct):    est=%.4f  95CI [%.4f, %.4f]  p=%.3f\n",
              s$z0, s$z0.ci[1], s$z0.ci[2], s$z0.p))
  cat(sprintf("  Total effect:     est=%.4f  95CI [%.4f, %.4f]  p=%.3f\n",
              s$tau.coef, s$tau.ci[1], s$tau.ci[2], s$tau.p))
  cat(sprintf("  Prop. mediated:   est=%.3f  95CI [%.3f, %.3f]\n",
              s$n0, s$n0.ci[1], s$n0.ci[2]))
}

# ============================================================
# 3. M0 — Primary Mediation: Change → Arousal → Attraction
#    No covariates; random intercepts only
# ============================================================

## 3a. Path a: Change → Arousal (mediator model)
m0_a <- lme4::lmer(Arousal ~ Change + Change2 + (1 | id),
             data = s3, REML = FALSE)

## 3b. Path b/c': Arousal + Change → Attraction (outcome model)
m0_b <- lme4::lmer(Attraction ~ Arousal + Change + Change2 + (1 | id),
             data = s3, REML = FALSE)

cat("\n=== M0: Path a — Change → Arousal ===\n")
print(summary(m0_a))
cat(sprintf("  R²m=%.3f  R²c=%.3f\n",
            MuMIn::r.squaredGLMM(m0_a)[1], MuMIn::r.squaredGLMM(m0_a)[2]))

cat("\n=== M0: Path b — Arousal → Attraction (controlling Change) ===\n")
print(summary(m0_b))
cat(sprintf("  R²m=%.3f  R²c=%.3f\n",
            MuMIn::r.squaredGLMM(m0_b)[1], MuMIn::r.squaredGLMM(m0_b)[2]))

## 3c. Bootstrap indirect effect
set.seed(42)
med_m0 <- mediation::mediate(
  model.m       = m0_a,
  model.y       = m0_b,
  treat         = "Change",
  mediator      = "Arousal",
  treat.value   = treat_val,
  control.value = ctrl_val,
  boot          = FALSE,
  # boot.ci.type not used with lmer (quasi-Bayesian approximation)
  sims          = 1000
)
print_med(med_m0, "M0: Primary mediation (no covariates)")

# ============================================================
# 4. M1 — Moderated Mediation: Accuracy moderates path a
#
# Theoretical prediction: Change → Arousal path is stronger on
# accurate (detected) trials — awareness gates misattribution.
# Path b (Arousal → Attraction) should not differ by Accuracy.
#
# Two approaches:
#   4a. Subset: separate mediation models for hits and misses
#       (simpler, but note selection bias: larger |Change| → higher Accuracy)
#   4b. Interaction: Change × Accuracy on path a; conditional ACME
# ============================================================

## ── 4a. Subset approach ──────────────────────────────────────

s3_hits   <- dplyr::filter(s3, Accuracy == 1)
s3_misses <- dplyr::filter(s3, Accuracy == 0)

# Path models for hits
m1_a_hit <- lme4::lmer(Arousal ~ Change + Change2 + (1 | id),
                 data = s3_hits, REML = FALSE)
m1_b_hit <- lme4::lmer(Attraction ~ Arousal + Change + Change2 + (1 | id),
                 data = s3_hits, REML = FALSE)

# Path models for misses
m1_a_miss <- lme4::lmer(Arousal ~ Change + Change2 + (1 | id),
                  data = s3_misses, REML = FALSE)
m1_b_miss <- lme4::lmer(Attraction ~ Arousal + Change + Change2 + (1 | id),
                  data = s3_misses, REML = FALSE)

# Use same contrast values but capped to observed data range per subset
treat_hit  <- min(treat_val,  max(s3_hits$Change))
ctrl_hit   <- max(ctrl_val,   min(s3_hits$Change))
treat_miss <- min(treat_val,  max(s3_misses$Change))
ctrl_miss  <- max(ctrl_val,   min(s3_misses$Change))

set.seed(42)
med_hits <- mediation::mediate(
  model.m = m1_a_hit, model.y = m1_b_hit,
  treat = "Change", mediator = "Arousal",
  treat.value = treat_hit, control.value = ctrl_hit,
  boot = FALSE, sims = 1000
)

set.seed(42)
med_misses <- mediation::mediate(
  model.m = m1_a_miss, model.y = m1_b_miss,
  treat = "Change", mediator = "Arousal",
  treat.value = treat_miss, control.value = ctrl_miss,
  boot = FALSE, sims = 1000
)

print_med(med_hits,   "M1a: Mediation on ACCURATE trials (hits)")
print_med(med_misses, "M1a: Mediation on INACCURATE trials (misses)")

cat("\n[Note] Subset approach caveat: larger |Change| increases both Accuracy",
    "and Arousal. Difference in ACME between hits and misses partly reflects",
    "change-magnitude differences, not pure awareness gating.",
    "See M1b interaction model for the principled test.\n")

## ── 4b. Interaction approach ─────────────────────────────────
# Change × Accuracy on path a; Accuracy (no interaction) on path b.
# Accuracy treated as continuous (0/1) for model compatibility.

m1_a_int <- lmerTest::lmer(Arousal ~ Change * Accuracy + Change2 + (1 | id),
                 data = s3, REML = FALSE)
m1_b_int <- lme4::lmer(Attraction ~ Arousal + Change + Change2 + Accuracy + (1 | id),
                 data = s3, REML = FALSE)

cat("\n=== M1b: Path a — Change × Accuracy interaction on Arousal ===\n")
print(summary(m1_a_int))
cat(sprintf("  R²m=%.3f  R²c=%.3f\n",
            MuMIn::r.squaredGLMM(m1_a_int)[1], MuMIn::r.squaredGLMM(m1_a_int)[2]))

# LRT: does adding Accuracy × Change improve path a over main effects?
m1_a_main <- lmerTest::lmer(Arousal ~ Change + Accuracy + Change2 + (1 | id),
                  data = s3, REML = FALSE)
cat("\n=== LRT: Accuracy × Change interaction on path a ===\n")
print(anova(m1_a_main, m1_a_int))

# Conditional indirect effects computed analytically from m1_a_int and m1_b_int.
# Refitting with fixed Accuracy causes rank deficiency (Change*constant = Change),
# so we extract coefficients directly instead.
#
# ACME at Accuracy=k = (b_Change + b_Change:Accuracy * k) * b_Arousal * (treat - ctrl)

cf_a <- fixef(m1_a_int)
cf_b <- fixef(m1_b_int)

b_change      <- cf_a["Change"]
b_interaction <- cf_a["Change:Accuracy"]
b_arousal     <- cf_b["Arousal"]
window        <- treat_val - ctrl_val

acme_miss <- b_change * b_arousal * window           # Accuracy = 0
acme_hit  <- (b_change + b_interaction) * b_arousal * window  # Accuracy = 1

cat(sprintf(
  "\nM1b: Conditional indirect effects (analytic):\n"
))
cat(sprintf("  Accuracy=0 (misses): ACME = %.4f\n", acme_miss))
cat(sprintf("  Accuracy=1 (hits):   ACME = %.4f\n", acme_hit))
cat(sprintf("  Difference (hits - misses): %.4f\n", acme_hit - acme_miss))
cat(sprintf(
  "\n  Note: SEs/CIs for conditional ACMEs from Bayesian M1 model (Section 7b)\n"
))

# ============================================================
# 5. M2 — Fatigue Sensitivity
#    Add BlockNum (trial within run) and Run_num to both path models
# ============================================================
cat("\n\n=== M2: Fatigue sensitivity — adding BlockNum + Run ===\n")

m2_a <- lme4::lmer(Arousal ~ Change + Change2 + BlockNum + Run_num + (1 | id),
             data = s3, REML = FALSE)
m2_b <- lme4::lmer(Attraction ~ Arousal + Change + Change2 + BlockNum + Run_num + (1 | id),
             data = s3, REML = FALSE)

cat("\nFatigue path a:\n"); print(summary(m2_a))
cat("\nFatigue path b:\n"); print(summary(m2_b))

set.seed(42)
med_m2 <- mediation::mediate(
  model.m = m2_a, model.y = m2_b,
  treat = "Change", mediator = "Arousal",
  treat.value = treat_val, control.value = ctrl_val,
  boot = FALSE, sims = 1000
)
print_med(med_m2, "M2: Primary mediation + fatigue covariates")

# Compare ACME to M0 — does controlling for fatigue meaningfully shrink it?
cat("\nM0 ACME:", round(summary(med_m0)$d0, 4),
    "| M2 ACME:", round(summary(med_m2)$d0, 4),
    "| Difference:", round(summary(med_m2)$d0 - summary(med_m0)$d0, 4), "\n")

# ============================================================
# 6. M3 — Salience Sensitivity
#    Add Salience to both path models (over and above fatigue)
# ============================================================
cat("\n\n=== M3: Salience sensitivity ===\n")

# Salience already canonical ("High"/"Low"); treatment-code with Low as reference
s3 <- s3 |>
  dplyr::mutate(Salience = factor(Salience, levels = c("Low", "High")))

m3_a <- lme4::lmer(Arousal ~ Change + Change2 + Salience + BlockNum + Run_num + (1 | id),
             data = s3, REML = FALSE)
m3_b <- lme4::lmer(Attraction ~ Arousal + Change + Change2 + Salience + BlockNum + Run_num + (1 | id),
             data = s3, REML = FALSE)

cat("\nSalience path a:\n"); print(summary(m3_a))
cat("\nSalience path b:\n"); print(summary(m3_b))

set.seed(42)
med_m3 <- mediation::mediate(
  model.m = m3_a, model.y = m3_b,
  treat = "Change", mediator = "Arousal",
  treat.value = treat_val, control.value = ctrl_val,
  boot = FALSE, sims = 1000
)
print_med(med_m3, "M3: Mediation + fatigue + Salience")

# ============================================================
# 7. Bayesian Mediation (brms)
#    Multivariate formulation: simultaneous estimation of path a and path b
#    Correct multilevel bootstrap via posterior sampling
# ============================================================
cache_m0  <- file.path(cache_dir, "brms_m0.rds")
cache_m1  <- file.path(cache_dir, "brms_m1.rds")

brms_ctrl <- list(adapt_delta = 0.95, max_treedepth = 12)

## ── 7a. Bayesian M0: primary mediation ───────────────────────

bf_a0 <- brms::bf(Arousal    ~ Change + Change2 + (1 | id))
bf_b0 <- brms::bf(Attraction ~ Arousal + Change + Change2 + (1 | id))

if (file.exists(cache_m0)) {
  brm_m0 <- readRDS(cache_m0)
  message("Loaded cached brms M0")
} else {
  brm_m0 <- brms::brm(
    formula   = bf_a0 + bf_b0 + set_rescor(FALSE),
    data      = s3,
    chains    = 4, iter = 8000, warmup = 2000,
    cores     = 4, seed  = 42,
    control   = brms_ctrl,
    prior     = c(brms::prior(exponential(1), class = "sd", resp = "Arousal"),
                  brms::prior(exponential(1), class = "sd", resp = "Attraction")),
    file      = cache_m0
  )
}

## Compute posterior distribution of indirect effect
## ACME = (b_Change_on_Arousal) × (b_Arousal_on_Attraction)
## Scaled to the same 2-SD window as the frequentist models
post_m0 <- brms::as_draws_df(brm_m0)

post_m0 <- post_m0 |>
  dplyr::mutate(
    # Path a coefficient (Change → Arousal)
    a_change  = b_Arousal_Change,
    # Path b coefficient (Arousal → Attraction)
    b_arousal = b_Attraction_Arousal,
    # Indirect effect scaled to 2-SD window of Change
    acme_bayes = a_change * b_arousal * (treat_val - ctrl_val)
  )

acme_summary <- post_m0 |>
  dplyr::summarise(
    mean   = mean(acme_bayes),
    median = median(acme_bayes),
    lower  = quantile(acme_bayes, 0.025),
    upper  = quantile(acme_bayes, 0.975),
    p_pos  = mean(acme_bayes > 0)   # posterior probability of positive indirect effect
  )

cat("\n=== Bayesian M0: Indirect effect (ACME) ===\n")
print(acme_summary)

## ── 7b. Bayesian M1: moderated mediation (Accuracy × Change on path a) ──

bf_a1 <- brms::bf(Arousal    ~ Change * Accuracy + Change2 + (1 | id))
bf_b1 <- brms::bf(Attraction ~ Arousal + Change + Change2 + Accuracy + (1 | id))

if (file.exists(cache_m1)) {
  brm_m1 <- readRDS(cache_m1)
  message("Loaded cached brms M1")
} else {
  brm_m1 <- brms::brm(
    formula   = bf_a1 + bf_b1 + set_rescor(FALSE),
    data      = s3,
    chains    = 4, iter = 8000, warmup = 2000,
    cores     = 4, seed  = 42,
    control   = brms_ctrl,
    prior     = c(brms::prior(exponential(1), class = "sd", resp = "Arousal"),
                  brms::prior(exponential(1), class = "sd", resp = "Attraction")),
    file      = cache_m1
  )
}

## Conditional indirect effects at Accuracy = 0 and Accuracy = 1
## For moderated path a: a(Acc) = b_Change + b_Change:Accuracy × Acc
post_m1 <- brms::as_draws_df(brm_m1)

post_m1 <- post_m1 |>
  dplyr::mutate(
    b_arousal       = b_Attraction_Arousal,
    # Path a at Accuracy = 1 (hits)
    a_hits          = b_Arousal_Change + `b_Arousal_Change:Accuracy`,
    # Path a at Accuracy = 0 (misses)
    a_misses        = b_Arousal_Change,
    # Conditional ACME (scaled to 2-SD window)
    acme_hits       = a_hits   * b_arousal * (treat_val - ctrl_val),
    acme_misses     = a_misses * b_arousal * (treat_val - ctrl_val),
    acme_difference = acme_hits - acme_misses   # key contrast
  )

cond_summary <- post_m1 |>
  dplyr::summarise(
    acme_hits_mean   = mean(acme_hits),
    acme_hits_lower  = quantile(acme_hits, 0.025),
    acme_hits_upper  = quantile(acme_hits, 0.975),
    acme_hits_p_pos  = mean(acme_hits > 0),
    acme_miss_mean   = mean(acme_misses),
    acme_miss_lower  = quantile(acme_misses, 0.025),
    acme_miss_upper  = quantile(acme_misses, 0.975),
    acme_miss_p_pos  = mean(acme_misses > 0),
    diff_mean        = mean(acme_difference),
    diff_lower       = quantile(acme_difference, 0.025),
    diff_upper       = quantile(acme_difference, 0.975),
    diff_p_hits_gt   = mean(acme_difference > 0)   # P(ACME_hits > ACME_misses)
  )

cat("\n=== Bayesian M1: Conditional indirect effects ===\n")
print(t(cond_summary))

# ============================================================
# 8. Summary Table
# ============================================================
mediation_table <- tibble::tribble(
  ~Model,   ~Covariates,               ~ACME,                 ~ACME_lower,          ~ACME_upper,          ~p_ACME,
  "M0",     "None",                    summary(med_m0)$d0,    summary(med_m0)$d0.ci[1], summary(med_m0)$d0.ci[2], summary(med_m0)$d0.p,
  "M2",     "Fatigue",                 summary(med_m2)$d0,    summary(med_m2)$d0.ci[1], summary(med_m2)$d0.ci[2], summary(med_m2)$d0.p,
  "M3",     "Fatigue + Salience",      summary(med_m3)$d0,    summary(med_m3)$d0.ci[1], summary(med_m3)$d0.ci[2], summary(med_m3)$d0.p
)

cat("\n=== ACME across model specifications (robustness) ===\n")
print(mediation_table)

cat("\n=== M1a: Awareness gating (subset approach) ===\n")
cat(sprintf("  Hits:   ACME=%.4f [%.4f, %.4f]  p=%.3f\n",
            summary(med_hits)$d0, summary(med_hits)$d0.ci[1],
            summary(med_hits)$d0.ci[2], summary(med_hits)$d0.p))
cat(sprintf("  Misses: ACME=%.4f [%.4f, %.4f]  p=%.3f\n",
            summary(med_misses)$d0, summary(med_misses)$d0.ci[1],
            summary(med_misses)$d0.ci[2], summary(med_misses)$d0.p))

cat("\n=== Bayesian M1: Conditional ACME ===\n")
cat(sprintf("  Hits:       %.4f [%.4f, %.4f]  P(>0)=%.3f\n",
            cond_summary$acme_hits_mean, cond_summary$acme_hits_lower,
            cond_summary$acme_hits_upper, cond_summary$acme_hits_p_pos))
cat(sprintf("  Misses:     %.4f [%.4f, %.4f]  P(>0)=%.3f\n",
            cond_summary$acme_miss_mean, cond_summary$acme_miss_lower,
            cond_summary$acme_miss_upper, cond_summary$acme_miss_p_pos))
cat(sprintf("  Difference: %.4f [%.4f, %.4f]  P(hits>misses)=%.3f\n",
            cond_summary$diff_mean, cond_summary$diff_lower,
            cond_summary$diff_upper, cond_summary$diff_p_hits_gt))

# ============================================================
# 9. Figure — Mediation diagram for hits vs misses
# ============================================================
# Forest-style plot showing ACME and 95% CI for each model/condition

plot_data <- tibble::tribble(
  ~Condition,    ~Model,          ~est,                   ~lower,                     ~upper,
  "Hits",        "Frequentist",   summary(med_hits)$d0,   summary(med_hits)$d0.ci[1], summary(med_hits)$d0.ci[2],
  "Misses",      "Frequentist",   summary(med_misses)$d0, summary(med_misses)$d0.ci[1], summary(med_misses)$d0.ci[2],
  "Overall",     "Frequentist",   summary(med_m0)$d0,     summary(med_m0)$d0.ci[1],   summary(med_m0)$d0.ci[2],
  "Hits",        "Bayesian",      cond_summary$acme_hits_mean,  cond_summary$acme_hits_lower,  cond_summary$acme_hits_upper,
  "Misses",      "Bayesian",      cond_summary$acme_miss_mean,  cond_summary$acme_miss_lower,  cond_summary$acme_miss_upper,
  "Overall",     "Bayesian",      acme_summary$mean,      acme_summary$lower,          acme_summary$upper
) |>
  dplyr::mutate(
    Condition = factor(Condition, levels = c("Overall", "Hits", "Misses")),
    Model     = factor(Model, levels = c("Frequentist", "Bayesian"))
  )

p_forest <- ggplot(plot_data, aes(x = est, y = Condition, colour = Model,
                                   xmin = lower, xmax = upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(height = 0.2), position = position_dodge(width = 0.5),
                 linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  scale_colour_bcat() +
  labs(
    x      = "Indirect effect (ACME)\nChange → Arousal → Attraction",
    y      = NULL,
    colour = NULL,
    title  = "Mediation: Change → Arousal → Attraction",
    subtitle = sprintf("ACME contrasts Change at +1 SD vs −1 SD (2-SD window = %.2f units)",
                       treat_val - ctrl_val)
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "mediation_forest.pdf"), p_forest,
       width = 7, height = 4, device = "pdf")

cat("\nFigure saved:", file.path(fig_dir, "mediation_forest.pdf"), "\n")
cat("\nStudy 3 mediation analysis complete.\n")

# ============================================================
# 10. MAIA AS PROTECTION AGAINST MISATTRIBUTION
#
# High MAIA → stronger awareness → weaker misattribution effect?
# Tests whether MAIA_total_z (participant-level) moderates path a
# (Change → Arousal) using MAIA as a Level-2 predictor.
#
# Prediction: negative Change × MAIA interaction on Arousal —
# higher MAIA participants show weaker arousal response to change,
# because better metacognitive calibration allows them to identify
# the bodily source rather than attributing it elsewhere.
#
# Also tests whether MAIA's effect on misattribution is mediated
# by the within-person Awareness index: MAIA → Awareness →
# reduced arousal transfer (mediation of moderation).
# ============================================================
cat("\n\n========================================\n")
cat("10. MAIA AS PROTECTION AGAINST MISATTRIBUTION\n")
cat("(Exploratory — not preregistered)\n")
cat("========================================\n")

# Join MAIA_total_z from summary to trial data
s3l_maia <- s3 |>
  dplyr::left_join(
    s3s |>
      dplyr::select(id, MAIA_total, Awareness) |>
      dplyr::mutate(
        id           = as.integer(id),
        MAIA_total_z = scale(MAIA_total)[, 1]
      ),
    by = "id"
  ) |>
  dplyr::filter(!is.na(MAIA_total_z))

cat(sprintf("\nN with MAIA data: %d participants, %d trials\n",
            dplyr::n_distinct(s3l_maia$id), nrow(s3l_maia)))

# ── 10a. MAIA moderates path a (Change → Arousal) ─────────────
m10_a_base <- lmerTest::lmer(Arousal ~ Change + Change2          + (1 | id),
                              data = s3l_maia, REML = FALSE)
m10_a_maia <- lmerTest::lmer(Arousal ~ Change * MAIA_total_z + Change2 + (1 | id),
                              data = s3l_maia, REML = FALSE)
m10_b      <- lmerTest::lmer(Attraction ~ Arousal + Change + Change2 + MAIA_total_z +
                               (1 | id), data = s3l_maia, REML = FALSE)

cat("\nLRT: does MAIA moderate Change → Arousal (path a)?\n")
print(anova(m10_a_base, m10_a_maia))
cat("\nPath a with MAIA moderation:\n")
print(broom.mixed::tidy(m10_a_maia, effects = "fixed", conf.int = TRUE), n = Inf)

# ── 10b. Conditional ACME at low/high MAIA (analytic) ────────
# Refitting with fixed MAIA_total_z causes rank deficiency.
# Compute conditional ACMEs directly from m10_a_maia coefficients:
#   ACME at MAIA=z = (b_Change + b_Change:MAIA * z) * b_Arousal * window

cf_10a <- lme4::fixef(m10_a_maia)
cf_10b <- lme4::fixef(m10_b)

b_change_10    <- cf_10a["Change"]
b_maia_int     <- cf_10a[grep("Change.*MAIA|MAIA.*Change",
                               names(cf_10a), value = TRUE)[1]]
b_arousal_10   <- cf_10b["Arousal"]
window_10      <- treat_val - ctrl_val
maia_lo        <- -1
maia_hi        <-  1

for (mval in c(maia_lo, maia_hi)) {
  label    <- if (mval < 0) "Low MAIA (-1 SD)" else "High MAIA (+1 SD)"
  acme_val <- (b_change_10 + b_maia_int * mval) * b_arousal_10 * window_10
  cat(sprintf("  %s: ACME = %.4f\n", label, acme_val))
}
cat(sprintf("  Difference (Low - High MAIA): %.4f\n",
            (b_change_10 + b_maia_int * maia_lo)  * b_arousal_10 * window_10 -
            (b_change_10 + b_maia_int * maia_hi) * b_arousal_10 * window_10))
cat("  Note: posterior CIs for conditional ACMEs from Bayesian model below\n")

# ── 10c. Does Awareness mediate the MAIA → protection pathway? ─
# MAIA → Awareness (metacognitive calibration) → attenuated path a
# Partial test: does controlling for Awareness reduce the MAIA × Change interaction?
if ("Awareness" %in% names(s3l_maia) && !all(is.na(s3l_maia$Awareness))) {
  m10_a_aw <- lmerTest::lmer(
    Arousal ~ Change * MAIA_total_z + Change * Awareness + Change2 + (1 | id),
    data = s3l_maia |> dplyr::filter(!is.na(Awareness)),
    REML = FALSE
  )
  cat("\nWith Awareness added to path a (MAIA × Change mediation test):\n")
  print(broom.mixed::tidy(m10_a_aw, effects = "fixed", conf.int = TRUE), n = Inf)
  cat("  If Change×MAIA attenuates with Awareness included → Awareness mediates MAIA protection\n")
}

# ── 10d. Bayesian MAIA moderation ─────────────────────────────
cache_m_maia <- file.path(cache_dir, "brms_maia_moderation.rds")
bf_a_maia <- brms::bf(Arousal    ~ Change * MAIA_total_z + Change2 + (1 | id))
bf_b_maia <- brms::bf(Attraction ~ Arousal + Change + Change2 + MAIA_total_z +
                          (1 | id))

if (file.exists(cache_m_maia)) {
  brm_maia <- readRDS(cache_m_maia)
  message("Loaded cached brms MAIA moderation model")
} else {
  brm_maia <- brms::brm(
    formula = bf_a_maia + bf_b_maia + set_rescor(FALSE),
    data    = s3l_maia,
    chains  = 4, iter = 8000, warmup = 2000,
    cores   = 4, seed  = 42,
    control = brms_ctrl,
    prior   = c(brms::prior(exponential(1), class = "sd", resp = "Arousal"),
                brms::prior(exponential(1), class = "sd", resp = "Attraction")),
    file    = cache_m_maia
  )
}

post_maia <- brms::as_draws_df(brm_maia)
# Posterior of the Change × MAIA interaction on path a
int_col <- grep("Change.*MAIA|MAIA.*Change",
                names(post_maia), value = TRUE)
if (length(int_col) > 0) {
  b_int <- post_maia[[int_col[1]]]
  cat(sprintf(
    "\nBayesian Change × MAIA on Arousal: b = %.4f [%.4f, %.4f]  P(b<0) = %.3f\n",
    mean(b_int), quantile(b_int, 0.025), quantile(b_int, 0.975),
    mean(b_int < 0)
  ))
  cat("  P(b<0) > 0.95 → high MAIA reliably attenuates arousal transfer\n")
}

# ============================================================
# 11. SAVE ALL RESULTS TO CSV
# ============================================================
cat("\nSaving mediation results to CSV...\n")

# ── Primary mediation table (M0, M2, M3 robustness) ──────────
readr::write_csv(mediation_table,
                 file.path(RESULTS_DIR, "study3_mediation_primary.csv"))

# ── M1a: Subset approach (hits vs misses) ─────────────────────
m1a_table <- tibble::tibble(
  model     = c("Hits (Accuracy=1)", "Misses (Accuracy=0)"),
  acme      = c(summary(med_hits)$d0,    summary(med_misses)$d0),
  acme_lower = c(summary(med_hits)$d0.ci[1],  summary(med_misses)$d0.ci[1]),
  acme_upper = c(summary(med_hits)$d0.ci[2],  summary(med_misses)$d0.ci[2]),
  p_acme    = c(summary(med_hits)$d0.p,  summary(med_misses)$d0.p),
  ade       = c(summary(med_hits)$z0,    summary(med_misses)$z0),
  p_ade     = c(summary(med_hits)$z0.p,  summary(med_misses)$z0.p),
  total     = c(summary(med_hits)$tau.coef, summary(med_misses)$tau.coef),
  prop_med  = c(summary(med_hits)$n0,    summary(med_misses)$n0),
  note      = "Selection bias caveat: hits have larger |Change| than misses"
)
readr::write_csv(m1a_table,
                 file.path(RESULTS_DIR, "study3_mediation_m1a_subset.csv"))

# ── M1b: Analytic conditional ACMEs ───────────────────────────
m1b_table <- tibble::tibble(
  model      = c("Analytic: Accuracy=0 (misses)",
                 "Analytic: Accuracy=1 (hits)",
                 "Difference (hits - misses)"),
  acme       = c(acme_miss, acme_hit, acme_hit - acme_miss),
  method     = "Analytic from lmer coefficients",
  note       = "CIs from Bayesian M1 model (see study3_mediation_bayesian.csv)"
)
readr::write_csv(m1b_table,
                 file.path(RESULTS_DIR, "study3_mediation_m1b_analytic.csv"))

# ── Bayesian M0: indirect effect ──────────────────────────────
readr::write_csv(
  acme_summary |> dplyr::mutate(model = "Bayesian M0"),
  file.path(RESULTS_DIR, "study3_mediation_bayesian_m0.csv")
)

# ── Bayesian M1: conditional ACMEs ────────────────────────────
bayes_m1_table <- tibble::tibble(
  condition  = c("Hits (Accuracy=1)", "Misses (Accuracy=0)", "Difference"),
  acme_mean  = c(cond_summary$acme_hits_mean,
                 cond_summary$acme_miss_mean,
                 cond_summary$diff_mean),
  acme_lower = c(cond_summary$acme_hits_lower,
                 cond_summary$acme_miss_lower,
                 cond_summary$diff_lower),
  acme_upper = c(cond_summary$acme_hits_upper,
                 cond_summary$acme_miss_upper,
                 cond_summary$diff_upper),
  p_direction = c(cond_summary$acme_hits_p_pos,
                  cond_summary$acme_miss_p_pos,
                  cond_summary$diff_p_hits_gt),
  interpretation = c(
    "P(ACME > 0) on detected trials",
    "P(ACME > 0) on missed trials",
    "P(hits ACME > misses ACME) — key awareness-gating test"
  )
)
readr::write_csv(bayes_m1_table,
                 file.path(RESULTS_DIR, "study3_mediation_bayesian_m1.csv"))

# ── MAIA protection: path a moderation ────────────────────────
maia_mod_cf <- broom.mixed::tidy(m10_a_maia, effects = "fixed",
                                  conf.int = TRUE)
readr::write_csv(maia_mod_cf,
                 file.path(RESULTS_DIR, "study3_maia_moderation_path_a.csv"))

# Conditional ACMEs at +/- 1 SD MAIA
maia_acme_table <- tibble::tibble(
  maia_level = c("Low MAIA (-1 SD)", "High MAIA (+1 SD)", "Difference"),
  maia_z     = c(maia_lo, maia_hi, NA),
  acme       = c(
    (b_change_10 + b_maia_int * maia_lo) * b_arousal_10 * window_10,
    (b_change_10 + b_maia_int * maia_hi) * b_arousal_10 * window_10,
    (b_change_10 + b_maia_int * maia_lo) * b_arousal_10 * window_10 -
      (b_change_10 + b_maia_int * maia_hi) * b_arousal_10 * window_10
  ),
  method     = "Analytic from lmer coefficients",
  note       = "Bayesian CIs from brms_maia_moderation.rds"
)
readr::write_csv(maia_acme_table,
                 file.path(RESULTS_DIR, "study3_maia_protection_acme.csv"))

# ── Bayesian MAIA moderation ───────────────────────────────────
if (exists("post_maia") && length(int_col) > 0) {
  bayes_maia_table <- tibble::tibble(
    parameter   = "Change x MAIA_z (path a)",
    mean        = mean(post_maia[[int_col[1]]]),
    lower_95    = quantile(post_maia[[int_col[1]]], 0.025),
    upper_95    = quantile(post_maia[[int_col[1]]], 0.975),
    p_negative  = mean(post_maia[[int_col[1]]] < 0),
    interpretation = "P(b<0)>0.95 → high MAIA reliably attenuates arousal transfer"
  )
  readr::write_csv(bayes_maia_table,
                   file.path(RESULTS_DIR, "study3_maia_bayesian_moderation.csv"))
}

message("\nStudy 3 mediation CSVs written to: ", RESULTS_DIR)
message("Files saved:")
message("  study3_mediation_primary.csv         — M0, M2, M3 ACME table")
message("  study3_mediation_m1a_subset.csv      — hits vs misses subset approach")
message("  study3_mediation_m1b_analytic.csv    — analytic conditional ACMEs")
message("  study3_mediation_bayesian_m0.csv     — Bayesian primary indirect effect")
message("  study3_mediation_bayesian_m1.csv     — Bayesian conditional ACMEs")
message("  study3_maia_moderation_path_a.csv    — MAIA x Change coefficients")
message("  study3_maia_protection_acme.csv      — conditional ACMEs at +/-1 SD MAIA")
if (exists("post_maia") && length(int_col) > 0)
  message("  study3_maia_bayesian_moderation.csv  — Bayesian MAIA moderation")
