# Test Block Accuracy Analysis
# BCAT Studies 4 and 5 — test block detection performance
#
# Examines how salience and direction affect detection accuracy
# in the calibrated test block (both conditions presented at each
# participant's individually estimated high-salience threshold).
#
# Analyses:
#   1. Descriptives: Pc by study / group / session / salience / direction
#   2. Trial-level GLMM: accuracy ~ salience * direction (+ group, session)
#   3. Person-level d' (3AFC, exact numerical formula) ~ salience + direction
#   4. Direction asymmetry: Faster vs Slower Pc within salience levels
#   5. Study 4 stratified (Breath vs Visual, between-group)
#   6. Study 5 session-specific (ses1_Breath / ses1_Visual / ses2_Breath)
#        + condition x session interaction
#
# NOTE on bias:
#   Full SDT criterion (c) requires the raw response (Faster/Slower/No-Change),
#   not just accuracy (0/1). Direction asymmetry (Pc_Faster vs Pc_Slower)
#   is used here as a proxy for directional response bias.
#
# Assumes BASE_DIR is set in environment (e.g. from MainAnalysis.R).
# Sources utils.R for compute_test_dprime_3afc() and dprime_from_pc().
#
# Outputs:
#   test_block_accuracy_descriptives.csv — Pc by ses_cond/salience/direction
#   test_block_dprime_summary.csv        — group-level d' (3AFC) by condition
#   test_block_dprime_by_group.csv       — d' with n_trials and Pc (for Table ST5)
#   test_block_direction_asymmetry.csv   — Faster vs Slower d' difference
#   test_block_accuracy_glmm.csv         — pooled GLMM fixed effects
#   test_block_dprime_lmm.csv            — LMM on person-level d'
#   test_block_study4_accuracy.csv       — Study 4 simple effects by group
#   test_block_study5_accuracy.csv       — Study 5 simple effects by session-condition

# Source utils.R for shared helpers (%||%, apply_bh_correction, etc.)
source(file.path(ANALYSIS_DIR, "utils.R"))

# d' for 3AFC (exact numerical integration, matches compute_test_dprime_3afc
# in utils.R -- redefined here because dprime_from_pc is a local closure there)
pc_from_dprime <- function(d) {
  stats::integrate(
    function(x) stats::dnorm(x - d) * stats::pnorm(x)^2,
    lower = -10, upper = 10
  )$value
}

dprime_from_pc <- function(pc) {
  if (pc <= 1/3) return(0)
  upper <- 1
  while (pc_from_dprime(upper) < pc) {
    upper <- upper + 1
    if (upper > 20) return(20)
  }
  stats::uniroot(
    function(d) pc_from_dprime(d) - pc,
    lower = 0, upper = upper
  )$root
}

# Optimizer control -- use for all lmer/glmer calls
glmer_ctrl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
lmer_ctrl  <- lme4::lmerControl( optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# ============================================================
# 1. LOAD AND HARMONIZE DATA
# ============================================================
s4_raw <- readr::read_csv(file.path(DATA_DIR, "study4_test.csv"))
s5_raw <- readr::read_csv(file.path(DATA_DIR, "study5_test.csv"))

s4 <- s4_raw |>
  dplyr::filter(!flag_zero_variance) |>
  dplyr::mutate(
    study     = "S4",
    salience  = dplyr::if_else(grepl("^high", Condition, ignore.case = TRUE), "High", "Low"),
    direction = dplyr::if_else(grepl("Acc$", Condition), "Faster", "Slower"),
    group     = Group,
    ses       = "ses1"
  ) |>
  dplyr::select(id, study, ses, group, salience, direction, Accuracy)

s5 <- s5_raw |>
  dplyr::mutate(
    study     = "S5",
    salience  = Salience,
    direction = Direction,
    group     = dplyr::if_else(Condition == "breath", "Breath", "Visual"),
    ses       = ses
  ) |>
  dplyr::select(id, study, ses, group, salience, direction, Accuracy)

test <- dplyr::bind_rows(s4, s5) |>
  dplyr::rename(accuracy = Accuracy) |>
  dplyr::mutate(
    id          = as.character(id),   # lme4 requires character/factor for (1|id)
    accuracy    = as.integer(accuracy),
    salience_f  = factor(salience,  levels = c("High", "Low")),
    direction_f = factor(direction, levels = c("Faster", "Slower")),
    group_f     = factor(group,     levels = c("Breath", "Visual")),
    study_f     = factor(study,     levels = c("S4", "S5")),
    ses_f       = factor(ses,       levels = c("ses1", "ses2")),
    salience_num = dplyr::if_else(salience_f == "High", 1L, 0L),
    dir_num      = dplyr::if_else(direction_f == "Faster", 1L, 0L)
  )

# Study 5: derive original group assignment (ses1 condition)
s5_assign <- dplyr::filter(test, study == "S5", ses == "ses1") |>
  dplyr::distinct(id, original_group = group)

test <- test |>
  dplyr::left_join(s5_assign, by = "id") |>
  dplyr::mutate(
    ses_cond = dplyr::case_when(
      study == "S4"                        ~ paste0("S4_", group),
      study == "S5" & ses == "ses1" & group == "Breath" ~ "S5_ses1_Breath",
      study == "S5" & ses == "ses1" & group == "Visual" ~ "S5_ses1_Visual",
      study == "S5" & ses == "ses2"                     ~ "S5_ses2_Breath"
    ),
    ses_cond_f = factor(ses_cond, levels = c(
      "S4_Breath", "S4_Visual",
      "S5_ses1_Breath", "S5_ses1_Visual", "S5_ses2_Breath"
    ))
  )

cat("Trial counts by ses_cond, salience, direction:\n")
test |> dplyr::count(ses_cond_f, salience_f, direction_f) |> print(n = 40)
cat("\nParticipant counts by ses_cond:\n")
test |> dplyr::distinct(id, ses_cond_f) |> dplyr::count(ses_cond_f) |> print()

# ============================================================
# 2. DESCRIPTIVES: Pc BY CONDITION
# ============================================================
# Group-level Pc (proportion correct) by all key factors
desc_pc <- test |>
  dplyr::group_by(ses_cond_f, salience_f, direction_f) |>
  dplyr::summarise(
    n_trials = dplyr::n(),
    n_pp     = dplyr::n_distinct(id),
    Pc       = round(mean(accuracy, na.rm = TRUE), 3),
    .groups  = "drop"
  )
cat("\n--- Descriptives: Pc by condition ---\n"); print(desc_pc, n = 40)

# Person-level Pc for d-prime computation
person_pc <- test |>
  dplyr::group_by(id, study, study_f, ses_cond_f, salience_f, direction_f,
                  group_f, ses_f, original_group) |>
  dplyr::summarise(
    n_trials = dplyr::n(),
    Pc       = mean(accuracy, na.rm = TRUE),
    .groups  = "drop"
  ) |>
  dplyr::filter(!is.na(Pc))

# ============================================================
# 3. TRIAL-LEVEL GLMM: ACCURACY ~ SALIENCE * DIRECTION
# ============================================================
# Pooled across studies; study_f, group_f, ses_f as covariates.
# Primary interest: salience_f (H5-adjacent), direction_f (bias check),
# and their interaction.

m_acc_main <- lme4::glmer(
  accuracy ~ salience_f * direction_f + group_f + study_f + ses_f +
    (1 | id),
  data   = test,
  family = binomial,
  control = glmer_ctrl
)
cat("\n--- Pooled GLMM: accuracy ~ salience * direction ---\n")
print(summary(m_acc_main)$coefficients)

# Additive model (no interaction) — baseline for main-effect LRTs
m_acc_additive <- lme4::glmer(
  accuracy ~ salience_f + direction_f + group_f + study_f + ses_f +
    (1 | id),
  data = test, family = binomial, control = glmer_ctrl
)

# LRT: is salience main effect significant? (additive vs. no salience)
m_acc_nosal <- lme4::glmer(
  accuracy ~ direction_f + group_f + study_f + ses_f + (1 | id),
  data = test, family = binomial, control = glmer_ctrl
)
lrt_sal <- anova(m_acc_nosal, m_acc_additive)
cat("\nLRT salience main effect:\n"); print(lrt_sal)

# LRT: is direction main effect significant? (additive vs. no direction)
m_acc_nodir <- lme4::glmer(
  accuracy ~ salience_f + group_f + study_f + ses_f + (1 | id),
  data = test, family = binomial, control = glmer_ctrl
)
lrt_dir <- anova(m_acc_nodir, m_acc_additive)
cat("\nLRT direction main effect:\n"); print(lrt_dir)

# LRT: is the salience × direction interaction significant?
lrt_int <- anova(m_acc_additive, m_acc_main)
cat("\nLRT salience × direction interaction:\n"); print(lrt_int)

# ============================================================
# 4. PERSON-LEVEL d' (3AFC) ~ SALIENCE AND DIRECTION
# ============================================================
# dprime_from_pc() from utils.R: exact numerical integration
# for 3-alternative forced choice (integral of Phi(x)^2 * phi(x-d'))

# Apply dprime to person-level Pc
person_dprime <- person_pc |>
  dplyr::mutate(
    # Floor/ceiling correction: avoid Pc <= 1/3 (chance) or Pc = 1
    Pc_corrected = pmax(pmin(Pc, 1 - 1 / (2 * n_trials)),
                        1/3 + 1 / (2 * n_trials)),
    dprime = purrr::map_dbl(Pc_corrected, dprime_from_pc)
  )

# Summary: mean d' by condition
dprime_summary <- person_dprime |>
  dplyr::group_by(ses_cond_f, salience_f, direction_f) |>
  dplyr::summarise(
    n_pp        = dplyr::n(),
    mean_dprime = round(mean(dprime, na.rm = TRUE), 3),
    sd_dprime   = round(sd(dprime,   na.rm = TRUE), 3),
    .groups     = "drop"
  )
cat("\n--- d' summary by condition ---\n"); print(dprime_summary, n = 40)

# LMM on person-level d': salience + direction + group + study + session
m_dprime_main <- lmerTest::lmer(
  dprime ~ salience_f * direction_f + group_f + study_f + ses_f +
    (1 | id),
  data = person_dprime, REML = FALSE, control = lmer_ctrl
)
cat("\n--- LMM: d' ~ salience * direction ---\n")
print(summary(m_dprime_main)$coefficients)

# Effect size: salience on d' (partial r)
sal_t  <- summary(m_dprime_main)$coefficients["salience_fLow", "t value"]
sal_df <- summary(m_dprime_main)$coefficients["salience_fLow", "df"]
cat(sprintf("\n  salience effect on d': partial r = %.3f\n",
            sal_t / sqrt(sal_t^2 + sal_df)))

# ============================================================
# 5. DIRECTION ASYMMETRY (BIAS PROXY)
# ============================================================
# Without raw response data, we use P(correct | Faster) vs P(correct | Slower)
# within each salience level as an index of directional detectability asymmetry.
# A Faster > Slower asymmetry suggests easier detection of accelerations.

dir_asym <- person_dprime |>
  dplyr::select(id, ses_cond_f, salience_f, direction_f, dprime) |>
  tidyr::pivot_wider(
    names_from  = direction_f,
    values_from = dprime
  ) |>
  dplyr::mutate(
    asym = Faster - Slower   # positive = Faster easier to detect
  ) |>
  dplyr::filter(!is.na(asym))

asym_summary <- dir_asym |>
  dplyr::group_by(ses_cond_f, salience_f) |>
  dplyr::summarise(
    n        = dplyr::n(),
    mean_asym = round(mean(asym), 3),
    sd_asym   = round(sd(asym),   3),
    t         = round(mean(asym) / (sd(asym) / sqrt(dplyr::n())), 3),
    p         = round(2 * pt(-abs(mean(asym) / (sd(asym) / sqrt(dplyr::n()))),
                              df = dplyr::n() - 1), 4),
    .groups   = "drop"
  )
cat("\n--- Direction asymmetry (d'_Faster - d'_Slower): ---\n")
print(asym_summary, n = 20)

# LMM test: is direction asymmetry consistent across groups / sessions?
m_asym <- lmerTest::lmer(
  asym ~ salience_f + ses_cond_f + (1 | id),
  data = dir_asym, REML = FALSE, control = lmer_ctrl
)
cat("\n--- LMM: direction asymmetry ~ salience + ses_cond ---\n")
print(summary(m_asym)$coefficients)

# ============================================================
# 6. STUDY 4 STRATIFIED: Breath vs Visual group
# ============================================================
test_s4 <- dplyr::filter(test, study == "S4")

m_s4_acc <- lme4::glmer(
  accuracy ~ salience_f * direction_f * group_f + (1 | id),
  data    = test_s4, family = binomial, control = glmer_ctrl
)
cat("\n--- Study 4: accuracy ~ salience * direction * group ---\n")
print(summary(m_s4_acc)$coefficients)

# Simple effects within each S4 group
s4_simple_acc <- purrr::map_dfr(c("Breath", "Visual"), function(grp) {
  df_g <- dplyr::filter(test_s4, group_f == grp)
  m <- lme4::glmer(
    accuracy ~ salience_f * direction_f + (1 | id),
    data = df_g, family = binomial, control = glmer_ctrl
  )
  broom.mixed::tidy(m, effects = "fixed") |>
    dplyr::mutate(study = "S4", group = grp,
                  dplyr::across(where(is.numeric), \(x) round(x, 4)))
})
cat("\n--- Study 4 simple effects by group ---\n"); print(s4_simple_acc)

# d' by ses_cond (consistent grouping for S4 and S5 -- ses_cond_f covers all rows)
# Join with desc_pc to include n_trials and Pc alongside d'
dprime_by_group <- person_dprime |>
  dplyr::group_by(ses_cond_f, salience_f, direction_f) |>
  dplyr::summarise(
    n           = dplyr::n(),
    mean_dprime = round(mean(dprime, na.rm = TRUE), 3),
    sd_dprime   = round(sd(dprime,   na.rm = TRUE), 3),
    .groups     = "drop"
  ) |>
  dplyr::left_join(
    desc_pc |> dplyr::select(ses_cond_f, salience_f, direction_f, n_trials, n_pp, Pc),
    by = c("ses_cond_f", "salience_f", "direction_f")
  )
cat("\n--- d' by ses_cond / salience / direction ---\n"); print(dprime_by_group)

# ============================================================
# 7. STUDY 5 SESSION-SPECIFIC: ses1_Breath / ses1_Visual / ses2_Breath
#    + condition x session interaction
# ============================================================
test_s5 <- dplyr::filter(test, study == "S5") |>
  dplyr::mutate(
    original_group_f = factor(original_group, levels = c("Breath", "Visual")),
    ses_f2           = factor(ses, levels = c("ses1", "ses2"))
  )

# Four-way: does gating pattern differ by condition AND session?
m_s5_cond_ses <- lme4::glmer(
  accuracy ~ salience_f * direction_f * original_group_f * ses_f2 + (1 | id),
  data    = test_s5, family = binomial, control = glmer_ctrl
)
cat("\n--- Study 5: condition x session interaction (accuracy) ---\n")
print(summary(m_s5_cond_ses)$coefficients)

# LRT: does condition x session improve fit over additive model?
m_s5_additive <- lme4::glmer(
  accuracy ~ salience_f * direction_f + original_group_f + ses_f2 + (1 | id),
  data = test_s5, family = binomial, control = glmer_ctrl
)
lrt_s5_ses <- anova(m_s5_additive, m_s5_cond_ses)
cat("\n--- LRT: condition x session interaction ---\n"); print(lrt_s5_ses)

# Simple effects within each ses_cond
s5_simple_acc <- purrr::map_dfr(
  c("S5_ses1_Breath", "S5_ses1_Visual", "S5_ses2_Breath"),
  function(sc) {
    df_sc <- dplyr::filter(test_s5, ses_cond == sc)
    m <- lme4::glmer(
      accuracy ~ salience_f * direction_f + (1 | id),
      data = df_sc, family = binomial, control = glmer_ctrl
    )
    broom.mixed::tidy(m, effects = "fixed") |>
      dplyr::mutate(study = "S5", ses_cond = sc,
                    dplyr::across(where(is.numeric), \(x) round(x, 4)))
  }
)
cat("\n--- Study 5 simple effects by session-condition ---\n"); print(s5_simple_acc)

# d' by session-condition
dprime_s5 <- person_dprime |>
  dplyr::filter(study == "S5") |>
  dplyr::group_by(ses_cond_f, salience_f, direction_f) |>
  dplyr::summarise(n = dplyr::n(), mean_dprime = round(mean(dprime), 3),
                   sd_dprime = round(sd(dprime), 3), .groups = "drop")
cat("\n--- Study 5 d' by session-condition / salience / direction ---\n"); print(dprime_s5)

# ============================================================
# 8. SAVE RESULTS
# ============================================================
# Pooled GLMM fixed effects
pooled_glmm <- broom.mixed::tidy(m_acc_main, effects = "fixed") |>
  dplyr::mutate(model = "accuracy_salience_direction",
                dplyr::across(where(is.numeric), \(x) round(x, 4)))

# LMM on d'
dprime_lmm <- broom.mixed::tidy(m_dprime_main, effects = "fixed") |>
  dplyr::mutate(model = "dprime_salience_direction",
                dplyr::across(where(is.numeric), \(x) round(x, 4)))

readr::write_csv(desc_pc,
  file.path(RESULTS_DIR, "test_block_accuracy_descriptives.csv"))
readr::write_csv(dprime_summary,
  file.path(RESULTS_DIR, "test_block_dprime_summary.csv"))
readr::write_csv(asym_summary,
  file.path(RESULTS_DIR, "test_block_direction_asymmetry.csv"))
readr::write_csv(pooled_glmm,
  file.path(RESULTS_DIR, "test_block_accuracy_glmm.csv"))
readr::write_csv(dprime_lmm,
  file.path(RESULTS_DIR, "test_block_dprime_lmm.csv"))
readr::write_csv(s4_simple_acc,
  file.path(RESULTS_DIR, "test_block_study4_accuracy.csv"))
readr::write_csv(s5_simple_acc,
  file.path(RESULTS_DIR, "test_block_study5_accuracy.csv"))
readr::write_csv(dplyr::bind_rows(dprime_by_group),
  file.path(RESULTS_DIR, "test_block_dprime_by_group.csv"))

cat("\nDone. Files written to:", RESULTS_DIR, "\n")
cat("  test_block_accuracy_descriptives.csv  -- Pc by ses_cond/salience/direction\n")
cat("  test_block_dprime_summary.csv         -- group-level d' (3AFC)\n")
cat("  test_block_direction_asymmetry.csv    -- Faster vs Slower d' difference\n")
cat("  test_block_accuracy_glmm.csv          -- pooled GLMM fixed effects\n")
cat("  test_block_dprime_lmm.csv             -- LMM on person-level d'\n")
cat("  test_block_study4_accuracy.csv        -- S4 simple effects by group\n")
cat("  test_block_study5_accuracy.csv        -- S5 simple effects by session-condition\n")
cat("  test_block_dprime_by_group.csv        -- d' by study/group/salience/direction\n")
