# ============================================================
# Study 2 Clean Script — BART Impulsivity Study
#
# Input files (place in base_dir/RawData/):
#   study2_taskdata.xls       — raw trial-level log (N=199)
#   study2_questionnairedata.xlsx   — Prolific demographic metadata
#
# Output files (written to data_dir):
#   study2_long.csv        — Circle1 trial-level rows
#   study2_summary.csv     — one row per participant
#   study2_exclusions.csv  — exclusion log for all 199 participants
#
# ============================================================

# Set Up ---------
packages <- c("readxl", "tidyverse", "psych")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
options(readr.show_col_types = FALSE)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

# ============================================================
# 1. Load Raw Data
# ============================================================
raw  <- read_excel(file.path(raw_dir, "study2_taskdata.xls"), col_types = "text")
demo <- read_excel(file.path(raw_dir, "study2_questionnairedata.xlsx"))

# Coerce numeric columns
num_cols <- c("Arousal", "Confidence", "DetectACC", "TotalRateChange",
              "TrackACC", "StepSize", "TrialNum", "EventNum",
              "reversalCount", "badcount", "goodcount",
              "BARTtotalpoints", "BARTtotalpumps",
              "BAS_Drive", "BAS_FunSeeking", "BAS_RewardResponsivenes", "BIS",
              "GA7_total", "PHQ15_total", "PHQ9_total", "MAIA_total",
              "REI_EA", "REI_EE", "REI_RA", "REI_RE",
              paste0("MAIA_", sprintf("%02d", 1:37)))
for (col in intersect(num_cols, names(raw))) {
  raw[[col]] <- as.numeric(raw[[col]])
}

# Coerce demo age
demo$age <- as.numeric(demo$age)

message("Loaded: ", nrow(raw), " rows, ", n_distinct(raw$UserId), " participants")

# ============================================================
# 2. Exclusion Screening
# ============================================================
all_ids   <- unique(raw$UserId)
bart_ids  <- unique(raw$UserId[raw$Block == "BART"])

circle1_var <- raw |>
  dplyr::filter(Block == "Circle1") |>
  dplyr::group_by(UserId) |>
  dplyr::summarise(
    arousal_sd  = sd(Arousal,    na.rm = TRUE),
    conf_sd     = sd(Confidence, na.rm = TRUE),
    n_c1_trials = n(),
    max_event   = max(as.numeric(EventNum), na.rm = TRUE),
    .groups = "drop"
  )

exclusion_log <- tibble(UserId = all_ids) |>
  dplyr::left_join(circle1_var, by = "UserId") |>
  dplyr::mutate(
    has_bart          = UserId %in% bart_ids,
    flag_no_bart      = !has_bart,
    flag_zero_arousal = has_bart & !is.na(arousal_sd) & arousal_sd == 0,
    flag_zero_conf    = has_bart & !is.na(conf_sd)    & conf_sd    == 0,
    excluded          = flag_no_bart | flag_zero_arousal | flag_zero_conf,
    exclusion_reason  = dplyr::case_when(
      flag_no_bart      ~ "task_noncompletion_no_BART",
      flag_zero_arousal ~ "zero_variance_arousal",
      flag_zero_conf    ~ "zero_variance_confidence",
      TRUE              ~ "included"
    )
  )

# Note: prior original hardcoded list of 14 IDs is a subset of the
# 32 caught by flag_no_bart; this criterion is more complete and principled.
prior_omit <- c("5eea192942776c07ee28e00c","589a46af57995d0001a8a03a",
               "58d43a081877ea000115a2f1","58dc5764bd5abb000185ee7b",
               "5ad9fe1b80491c000142e065","5be415810241050001ca2318",
               "5b75ac687dc4b00001af4e50","5d4104de2a3ce2001bb52d12",
               "5ddafb0f2e4a05a5ff7759c9","5e10b8aa67f17079770dc73f",
               "5e36b149a617e5637898ab7e","5e503fa37cbf4308d1131e52",
               "5ea1dd85692b7c0c73166152","5ecd204e8f0812042642a2a4")
exclusion_log <- exclusion_log |>
  dplyr::mutate(in_prior_omit_list = UserId %in% prior_omit)

n_excluded <- sum(exclusion_log$excluded)
n_final    <- sum(!exclusion_log$excluded)
message("Excluded: ", n_excluded, " | Final N: ", n_final)
message("  - Task non-completion (no BART): ",
        sum(exclusion_log$flag_no_bart))
message("  - Zero-variance arousal/conf (additional): ",
        sum(exclusion_log$flag_zero_arousal | exclusion_log$flag_zero_conf))
message("Note: prior original list had 14 IDs; ",
        sum(exclusion_log$in_prior_omit_list & exclusion_log$excluded),
        " are captured by the current criterion")

included_ids <- exclusion_log |> dplyr::filter(!excluded) |> dplyr::pull(UserId)

# ============================================================
# 3. Join Demographics
# ============================================================
# Merge Prolific demo onto raw (age, Sex, etc.)
raw <- raw |>
  dplyr::left_join(
    demo |> dplyr::select(UserId, age, Sex,
                          `Country of Birth`, `Current Country of Residence`,
                          `Employment Status`, `Nationality`),
    by = "UserId"
  )

# ============================================================
# 4. Prepare Long-Format Data (Circle1 trials only)
# ============================================================
long <- raw |>
  dplyr::filter(Block == "Circle1", UserId %in% included_ids) |>
  dplyr::mutate(
    study_id = 2L,
    NumericChangeType = dplyr::recode(ChangeType,
      "faster"   = -1L,
      "slower"   =  1L,
      "nochange" =  0L
    ),
    ChangeVector = TotalRateChange * NumericChangeType
  ) |>
  dplyr::select(
    study_id, UserId, TrialNum, EventNum,
    ChangeType, NumericChangeType, TotalRateChange, ChangeVector,
    DetectACC, Arousal, Confidence, TrackACC, StepSize, reversalCount,
    age, Sex
  )

message("Long-format rows: ", nrow(long))

# ============================================================
# 5. Compute Summary Variables
# ============================================================

## 5a. Threshold (mean TotalRateChange for TrialNum > 21)
thresh <- long |>
  dplyr::filter(TrialNum > 21) |>
  dplyr::group_by(UserId) |>
  dplyr::summarise(thresh_c1 = mean(TotalRateChange, na.rm = TRUE),
                   .groups = "drop")

## 5b. Signal Detection Theory
sdt <- long |>
  dplyr::mutate(
    signal_trial = ChangeType %in% c("faster", "slower"),
    hit  = signal_trial &  (DetectACC == 1),
    miss = signal_trial & !(DetectACC == 1),
    fa   = !signal_trial & !(DetectACC == 1),
    cr   = !signal_trial &  (DetectACC == 1)
  ) |>
  dplyr::group_by(UserId) |>
  dplyr::summarise(
    n_hits   = sum(hit,  na.rm = TRUE),
    n_misses = sum(miss, na.rm = TRUE),
    n_fa     = sum(fa,   na.rm = TRUE),
    n_cr     = sum(cr,   na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    H_corr  = (n_hits + 0.5) / (n_hits + n_misses + 1),
    FA_corr = (n_fa   + 0.5) / (n_fa   + n_cr     + 1),
    dprime  = qnorm(H_corr) - qnorm(FA_corr),
    c_bias  = -0.5 * (qnorm(H_corr) + qnorm(FA_corr))
  )

## 5c. Mean arousal, confidence, awareness
task_summary <- long |>
  dplyr::group_by(UserId) |>
  dplyr::summarise(
    mean_Arousal    = mean(Arousal,    na.rm = TRUE),
    mean_Confidence = mean(Confidence, na.rm = TRUE),
    Awareness       = suppressWarnings(cor(Confidence, DetectACC, use = "complete.obs")),
    .groups = "drop"
  )

## 5d. Questionnaire totals (from questionnaire_PHQ15 block which carries
##     all summary scores)
quest_cols <- c("BARTtotalpoints", "BARTtotalpumps",
                "BAS_Drive", "BAS_FunSeeking", "BAS_RewardResponsivenes", "BIS",
                "GA7_total", "PHQ15_total", "PHQ9_total", "MAIA_total",
                "REI_EA", "REI_EE", "REI_RA", "REI_RE",
                "MAIA_attentionregulation", "MAIA_bodylisten", "MAIA_emoaware",
                "MAIA_notdistracting", "MAIA_noticing", "MAIA_notworrying",
                "MAIA_selfreg", "MAIA_trusting")
quest <- raw |>
  dplyr::filter(Block == "questionnaire_PHQ15", UserId %in% included_ids) |>
  dplyr::group_by(UserId) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(UserId, all_of(intersect(quest_cols, names(raw))))

## 5e. Demographics
demo_summary <- raw |>
  dplyr::filter(UserId %in% included_ids) |>
  dplyr::group_by(UserId) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(UserId, Age, Gender, age, Sex)

## 5f. Join
summary_df <- demo_summary |>
  dplyr::left_join(thresh,       by = "UserId") |>
  dplyr::left_join(sdt,          by = "UserId") |>
  dplyr::left_join(task_summary, by = "UserId") |>
  dplyr::left_join(quest,        by = "UserId") |>
  dplyr::mutate(study_id = 2L) |>
  dplyr::select(study_id, UserId, everything())

message("Summary rows: ", nrow(summary_df), " (should be ", n_final, ")")

# ============================================================
# 5b. Standardise Output Column Names
# ============================================================
long_out <- long |>
  dplyr::rename(
    id       = UserId,
    Accuracy = DetectACC,
    Change   = ChangeVector
  ) |>
  dplyr::mutate(
    Direction = dplyr::case_when(
      ChangeType == "faster"   ~ "Faster",
      ChangeType == "slower"   ~ "Slower",
      ChangeType == "nochange" ~ "NoChange",
      TRUE ~ ChangeType
    )
  ) |>
  dplyr::select(-ChangeType)

summary_out <- summary_df |>
  dplyr::rename(
    id                  = UserId,
    MAIA_Noticing       = MAIA_noticing,
    MAIA_NotDistracting = MAIA_notdistracting,
    MAIA_NotWorrying    = MAIA_notworrying,
    MAIA_AttentionReg   = MAIA_attentionregulation,
    MAIA_EmoAware       = MAIA_emoaware,
    MAIA_SelfReg        = MAIA_selfreg,
    MAIA_BodyListen     = MAIA_bodylisten,
    MAIA_Trusting       = MAIA_trusting
  )

# ============================================================
# 6. Write Outputs
# ============================================================
write_csv(long_out,      file.path(data_dir, "study2_long.csv"))
write_csv(summary_out,   file.path(data_dir, "study2_summary.csv"))
write_csv(exclusion_log, file.path(data_dir, "study2_exclusions.csv"))

message("Study 2 cleaning complete.")
message("  Long:       ", nrow(long_out),    " rows")
message("  Summary:    ", nrow(summary_out), " participants")
message("  Exclusions: ", n_excluded, " removed (", n_final, " retained)")
message("  Note: Revised from prior N=169 → N=", n_final,
        " using principled non-completion criterion")
