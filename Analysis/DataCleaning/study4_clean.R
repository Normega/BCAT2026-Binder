# ============================================================
# Study 4 Clean Script — Online BCAT Preparation Study
# 
# Input files (place in base_dir/RawData/):
#   study4_longData.csv           — staircase trial-level data (N=131)
#   study4_testData.csv           — test-phase trial-level data (N=131)
#   study4_questionnairedata.xlsx — questionnaire data (N=~175 usable)
#
# Output files (written to data_dir):
#   study4_long.csv        — staircase trial-level rows
#   study4_test.csv        — test-phase trial-level rows
#   study4_summary.csv     — one row per participant
#   study4_exclusions.csv  — exclusion/flag log
#
# Exclusion policy:
#   Zero-variance Arousal or Confidence across staircase trials
#   are FLAGGED but NOT excluded (per pre-analysis decision).
#   Flag column `flag_zero_variance = TRUE` included throughout.
#   N retained = 131 (5 flagged).
#
# Questionnaire scoring:
#   MAIA-2 (24 items, 0–4): reverse items 4–9 (1-indexed)
#   SPANE:  positive items 1,3,5,7,10,12 | negative items 2,4,6,8,9,11
#   PWB:    mean of 8 items (1–7)
#   BIPS:   sum of 9 items (1–5); item 8 reverse-scored
#   BARQ:   mean of 12 items (0–3)
#   PHQ-4:  anxiety = items 1+2 | depression = items 3+4
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
staircase <- read_csv(file.path(raw_dir, "study4_longData.csv"))
test_phase <- read_csv(file.path(raw_dir, "study4_testData.csv"))
qualtrics  <- read_excel(file.path(raw_dir, "study4_questionnairedata.xlsx"),
                         sheet = "Sheet0") |>  # Sheet0 = questionnaire data;
                                                # Sheet1 is a participant tracking log
  # Standard Qualtrics export has 2 header rows:
  #   Row 0 (read as column names): machine-readable names e.g. SPANE_1, MAIA24_1
  #   Row 1 (first data row):       verbose question labels — drop this
  dplyr::slice(-1)

# Clean column names (trim whitespace, e.g. "Comprehension Check ")
names(qualtrics) <- trimws(names(qualtrics))

# Rename ambiguous Q columns immediately after loading.
# readxl disambiguates duplicate column names by appending ...N (position).
# The premood ratings (Q1/Q3/Q5/Q7) and age (Q2.1) appear as:
#   Q1       (pos 25)  — unique, kept as-is
#   Q3       (pos 27)  — unique, kept as-is
#   Q5...29  (pos 29)  — premood stress; Q5...111 is indigenous identity follow-up
#   Q7       (pos 31)  — unique, kept as-is
#   Q2...102 (pos 102) — age question (was Q2.1 in earlier file version)
qualtrics <- qualtrics |>
  dplyr::rename(
    premood_energy  = Q1,
    premood_valence = Q3,
    premood_stress  = `Q5...29`,
    premood_life    = Q7,
    age_raw         = `Q2...102`
  )

# Coerce SONAid to integer to match task data id
qualtrics <- qualtrics |>
  dplyr::mutate(id = suppressWarnings(as.integer(SONAid))) |>
  dplyr::filter(!is.na(id))

message("Staircase:  ", nrow(staircase),  " rows, ", n_distinct(staircase$id),  " participants")
message("Test phase: ", nrow(test_phase),  " rows, ", n_distinct(test_phase$id), " participants")
message("Qualtrics:  ", nrow(qualtrics),   " rows, ", n_distinct(qualtrics$id),  " participants")

# ============================================================
# 2. Score Questionnaires
# ============================================================

## MAIA-2 (24 items, scored 0–4)
maia24_cols    <- paste0("MAIA24_", 1:24)
maia24_reverse <- 4:9   # 1-indexed items to reverse (Not-Distracting: 4-6; Not-Worrying: 7-9)
maia24_max     <- 4L

# Qualtrics exported text labels rather than numeric values.
# Recode to 0–4 before scoring. These five labels are the only values
# present in the raw data (verified from study4_questionnairedata.xlsx).
maia_label_map <- c(
  "Never"               = 0L,
  "Sometimes"           = 1L,
  "About half the time" = 2L,
  "Most of the time"    = 3L,
  "Always"              = 4L
)

qualtrics <- qualtrics |>
  dplyr::mutate(across(
    all_of(maia24_cols),
    ~ dplyr::recode(.x, !!!maia_label_map, .default = NA_integer_)
  )) |>
  dplyr::mutate(across(all_of(maia24_cols), as.numeric)) |>
  dplyr::mutate(
    across(all_of(paste0("MAIA24_", maia24_reverse)),
           ~ maia24_max - .x)
  ) |>
  dplyr::mutate(
    # Sum scores (0–4 per item after reversal)
    MAIA_Noticing       = rowSums(across(paste0("MAIA24_", 1:3)),   na.rm = TRUE),  # max 12
    MAIA_NotDistracting = rowSums(across(paste0("MAIA24_", 4:6)),   na.rm = TRUE),  # max 12
    MAIA_NotWorrying    = rowSums(across(paste0("MAIA24_", 7:9)),   na.rm = TRUE),  # max 12
    MAIA_AttentionReg   = rowSums(across(paste0("MAIA24_", 10:13)), na.rm = TRUE),  # max 16
    MAIA_EmoAware       = rowSums(across(paste0("MAIA24_", 14:16)), na.rm = TRUE),  # max 12
    MAIA_SelfReg        = rowSums(across(paste0("MAIA24_", 17:18)), na.rm = TRUE),  # max  8
    MAIA_BodyListen     = rowSums(across(paste0("MAIA24_", 19:21)), na.rm = TRUE),  # max 12
    MAIA_Trusting       = rowSums(across(paste0("MAIA24_", 22:24)), na.rm = TRUE),  # max 12
    MAIA_total          = rowSums(across(all_of(maia24_cols)),      na.rm = TRUE)   # max 96
  )

## SPANE (12 items, 1–5)
# Positive affect items: 1,3,5,7,10,12
# Negative affect items: 2,4,6,8,9,11
# Qualtrics labels verified from raw file:
spane_label_map <- c(
  "Very rarely or never" = 1L,
  "Rarely"               = 2L,
  "Sometimes"            = 3L,
  "Often"                = 4L,
  "Very often or always" = 5L
)
spane_cols <- paste0("SPANE_", 1:12)
qualtrics <- qualtrics |>
  dplyr::mutate(across(
    all_of(spane_cols),
    ~ as.numeric(dplyr::recode(.x, !!!spane_label_map, .default = NA_integer_))
  )) |>
  dplyr::mutate(
    # Positive items 1,3,5,7,10,12 | Negative items 2,4,6,8,9,11; max = 30 each
    SPANE_Pos = rowSums(across(paste0("SPANE_", c(1,3,5,7,10,12))), na.rm = TRUE),
    SPANE_Neg = rowSums(across(paste0("SPANE_", c(2,4,6,8,9,11))),  na.rm = TRUE)
  )

## PWB (8 items, 1–7)
# Qualtrics labels verified from raw file (7-point Likert):
pwb_label_map <- c(
  "Strongly disagree"          = 1L,
  "Disagree"                   = 2L,
  "Somewhat disagree"          = 3L,
  "Neither agree nor disagree" = 4L,
  "Somewhat agree"             = 5L,
  "Agree"                      = 6L,
  "Strongly agree"             = 7L
)
pwb_cols <- paste0("PWB_", 1:8)
qualtrics <- qualtrics |>
  dplyr::mutate(across(
    all_of(pwb_cols),
    ~ as.numeric(dplyr::recode(.x, !!!pwb_label_map, .default = NA_integer_))
  )) |>
  dplyr::mutate(PWB_total = rowSums(across(all_of(pwb_cols)), na.rm = TRUE))  # max 56

## BIPS (9 items, 1–5; item 8 reverse-scored)
# Qualtrics labels verified from raw file:
bips_label_map <- c(
  "Never"        = 1L,
  "Almost never" = 2L,
  "Sometimes"    = 3L,
  "Fairly Often" = 4L,
  "Very Often"   = 5L
)
bips_cols <- paste0("BIPS_", 1:9)
qualtrics <- qualtrics |>
  dplyr::mutate(across(
    all_of(bips_cols),
    ~ as.numeric(dplyr::recode(.x, !!!bips_label_map, .default = NA_integer_))
  )) |>
  dplyr::mutate(BIPS_8 = 6L - BIPS_8) |>   # reverse item 8 (max 5 + 1 = 6)
  dplyr::mutate(BIPS_total = rowSums(across(all_of(bips_cols)), na.rm = TRUE))

## BARQ (12 items, 0–3)
# Qualtrics labels verified from raw file:
barq_label_map <- c(
  "Completely disagree" = 0L,
  "Somewhat disagree"   = 1L,
  "Somewhat agree"      = 2L,
  "Completely agree"    = 3L
)
barq_cols <- paste0("BARQ_", 1:12)
qualtrics <- qualtrics |>
  dplyr::mutate(across(
    all_of(barq_cols),
    ~ as.numeric(dplyr::recode(.x, !!!barq_label_map, .default = NA_integer_))
  )) |>
  dplyr::mutate(BARQ_total = rowSums(across(all_of(barq_cols)), na.rm = TRUE))  # max 36

## PHQ-4 (4 items, 0–3)
# Qualtrics labels verified from raw file:
phq4_label_map <- c(
  "Not at all"             = 0L,
  "Several days"           = 1L,
  "More than half the days" = 2L,
  "Nearly everyday"        = 3L
)
phq4_cols <- paste0("PHQ4_", 1:4)
qualtrics <- qualtrics |>
  dplyr::mutate(across(
    all_of(phq4_cols),
    ~ as.numeric(dplyr::recode(.x, !!!phq4_label_map, .default = NA_integer_))
  )) |>
  dplyr::mutate(
    PHQ4_Anxiety    = PHQ4_1 + PHQ4_2,
    PHQ4_Depression = PHQ4_3 + PHQ4_4
  )

## Pre-task mood ratings (columns renamed at load time to avoid NSE conflicts)
qualtrics <- qualtrics |>
  dplyr::mutate(
    premood_energy  = as.numeric(premood_energy),
    premood_valence = as.numeric(premood_valence),
    premood_stress  = as.numeric(premood_stress),
    premood_life    = as.numeric(premood_life),
    age             = as.numeric(age_raw)
  )

## Gender (Q3_1–Q3_6: check-all-that-apply)
# Qualtrics exports one column per option; each contains the label text if
# selected, NA if not. Collapse into a single ordered character variable,
# handling multi-selections (e.g., Non-binary + Woman) before single checks.
qualtrics <- qualtrics |>
  dplyr::mutate(
    Gender = dplyr::case_when(
      !is.na(Q3_3) & !is.na(Q3_2) ~ "Woman + Non-binary",
      !is.na(Q3_3) & !is.na(Q3_1) ~ "Man + Non-binary",
      !is.na(Q3_3)                 ~ "Non-binary",
      !is.na(Q3_4) & !is.na(Q3_2) ~ "Transgender Woman",
      !is.na(Q3_4) & !is.na(Q3_1) ~ "Transgender Man",
      !is.na(Q3_2)                 ~ "Woman",
      !is.na(Q3_1)                 ~ "Man",
      !is.na(Q3_5)                 ~ dplyr::coalesce(as.character(Q3_5_TEXT), "Not listed"),
      !is.na(Q3_6)                 ~ "Prefer not to answer",
      TRUE                         ~ NA_character_
    )
  )

message("Gender distribution (pre-exclusion): ",
        paste(names(table(qualtrics$Gender, useNA = "ifany")),
              table(qualtrics$Gender, useNA = "ifany"),
              sep = "=", collapse = ", "))

# ============================================================
# 3. Exclusion Screening
# ============================================================
all_ids <- unique(staircase$id)

arousal_var <- staircase |>
  dplyr::group_by(id) |>
  dplyr::summarise(
    arousal_sd = sd(as.numeric(Arousal),    na.rm = TRUE),
    conf_sd    = sd(as.numeric(Confidence), na.rm = TRUE),
    n_trials   = n(),
    .groups = "drop"
  )

exclusion_log <- tibble(id = all_ids) |>
  dplyr::left_join(arousal_var, by = "id") |>
  dplyr::mutate(
    has_questionnaire  = id %in% qualtrics$id,
    flag_zero_arousal  = !is.na(arousal_sd) & arousal_sd == 0,
    flag_zero_conf     = !is.na(conf_sd)    & conf_sd    == 0,
    flag_zero_variance = flag_zero_arousal | flag_zero_conf,
    excluded           = FALSE,    # flagged only, not excluded
    exclusion_reason   = dplyr::case_when(
      flag_zero_variance ~ "zero_variance_flagged_retained",
      !has_questionnaire ~ "no_questionnaire_match",
      TRUE               ~ "included"
    )
  )

message("Total N: ", length(all_ids))
message("Zero-variance flagged (retained): ", sum(exclusion_log$flag_zero_variance))
message("No questionnaire match: ", sum(!exclusion_log$has_questionnaire))

# ============================================================
# 4. Prepare Long-Format Staircase Data
# ============================================================
long <- staircase |>
  dplyr::left_join(
    exclusion_log |> dplyr::select(id, flag_zero_variance, exclusion_reason),
    by = "id"
  ) |>
  dplyr::mutate(
    study_id  = 4L,
    Level     = as.numeric(Level),
    Accuracy  = as.numeric(Accuracy),
    Arousal   = as.numeric(Arousal),
    Confidence = as.numeric(Confidence),
    Change    = as.numeric(Change),
    Salience  = factor(Salience, levels = c("Low", "High")),
    Direction = factor(Direction, levels = c("Faster", "NoChange", "Slower")),
    Group     = factor(Group, levels = c("Breath", "Visual"))
  ) |>
  dplyr::select(
    study_id, id, Group, Trial, Condition, Salience, Direction, TrialDirection,
    Level, Change, Accuracy, Confidence, Arousal,
    flag_zero_variance, exclusion_reason
  )

# Prepare test-phase data
test_out <- test_phase |>
  dplyr::left_join(
    exclusion_log |> dplyr::select(id, flag_zero_variance),
    by = "id"
  ) |>
  dplyr::mutate(
    study_id   = 4L,
    Level      = as.numeric(Level),
    Accuracy   = as.numeric(Accuracy),
    Arousal    = as.numeric(Arousal),
    Confidence = as.numeric(Confidence),
    # Recode numeric codes to string labels before factoring
    # (study4_testData.csv stores Salience as 0/1 and Direction as -1/1)
    Salience  = dplyr::case_when(
      Salience == 0 ~ "Low",
      Salience == 1 ~ "High",
      TRUE          ~ as.character(Salience)
    ),
    Direction = dplyr::case_when(
      Direction == -1 ~ "Faster",
      Direction ==  1 ~ "Slower",
      TRUE            ~ as.character(Direction)
    ),
    Salience   = factor(Salience, levels = c("Low", "High")),
    Direction  = factor(Direction, levels = c("Faster", "Slower")),
    Group      = factor(Group, levels = c("Breath", "Visual"))
  ) |>
  dplyr::group_by(id) |>
  dplyr::mutate(Trial = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::select(
    study_id, id, Group, Trial, Condition, Salience, Direction,
    Level, Accuracy, Confidence, Arousal, flag_zero_variance
  )

message("Staircase long rows: ", nrow(long))
message("Test phase rows:     ", nrow(test_out))

# ============================================================
# 5. Compute Summary Variables
# ============================================================

## 5a. Threshold per condition (mean Level for Trial >= 8)
thresh_staircase <- long |>
  dplyr::filter(as.numeric(Trial) >= 8) |>
  dplyr::group_by(id, Salience, Direction) |>
  dplyr::summarise(threshold = mean(Level, na.rm = TRUE), .groups = "drop") |>
  dplyr::filter(Direction != "NoChange") |>
  tidyr::pivot_wider(
    names_from  = c(Salience, Direction),
    values_from = threshold,
    names_glue  = "thresh_{Salience}_{Direction}"
  )

## 5b. Test-phase threshold (mean Level per condition)
thresh_test <- test_out |>
  dplyr::group_by(id, Salience, Direction) |>
  dplyr::summarise(test_threshold = mean(Level, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(
    names_from  = c(Salience, Direction),
    values_from = test_threshold,
    names_glue  = "test_thresh_{Salience}_{Direction}"
  )

## 5c. Signal detection (test phase: high-salience threshold for both conditions)
sdt <- long |>
  dplyr::filter(Direction != "NoChange") |>
  dplyr::mutate(
    signal_trial = TRUE,
    hit  = Accuracy == 1,
    miss = Accuracy != 1
  ) |>
  dplyr::bind_rows(
    long |>
      dplyr::filter(Direction == "NoChange") |>
      dplyr::mutate(signal_trial = FALSE,
                    fa = Accuracy != 1,
                    cr = Accuracy == 1)
  ) |>
  dplyr::group_by(id) |>
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

## 5d. Mean Arousal, Confidence, Awareness
task_summary <- long |>
  dplyr::group_by(id) |>
  dplyr::summarise(
    Group           = first(Group),
    mean_Arousal    = mean(Arousal,    na.rm = TRUE),
    mean_Confidence = mean(Confidence, na.rm = TRUE),
    Awareness       = suppressWarnings(cor(Confidence, Accuracy, use = "complete.obs")),
    n_trials        = n(),
    .groups = "drop"
  )

## 5e. Questionnaire summary
quest_summary <- qualtrics |>
  dplyr::filter(id %in% all_ids) |>
  dplyr::group_by(id) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(id, age, Condition, Gender,
                MAIA_total, MAIA_Noticing, MAIA_NotDistracting, MAIA_NotWorrying,
                MAIA_AttentionReg, MAIA_EmoAware, MAIA_SelfReg,
                MAIA_BodyListen, MAIA_Trusting,
                SPANE_Pos, SPANE_Neg, PWB_total, BIPS_total, BARQ_total,
                PHQ4_Anxiety, PHQ4_Depression,
                premood_energy, premood_valence, premood_stress, premood_life)

## 5f. Join
summary_df <- task_summary |>
  dplyr::left_join(thresh_staircase, by = "id") |>
  dplyr::left_join(thresh_test,      by = "id") |>
  dplyr::left_join(sdt,              by = "id") |>
  dplyr::left_join(quest_summary,    by = "id") |>
  dplyr::left_join(
    exclusion_log |> dplyr::select(id, flag_zero_variance, exclusion_reason),
    by = "id"
  ) |>
  dplyr::mutate(study_id = 4L) |>
  dplyr::select(study_id, id, Group, everything())

message("Summary rows: ", nrow(summary_df))

# ============================================================
# 6. Write Outputs
# ============================================================
write_csv(long,          file.path(data_dir, "study4_long.csv"))
write_csv(test_out,      file.path(data_dir, "study4_test.csv"))
write_csv(summary_df,    file.path(data_dir, "study4_summary.csv"))
write_csv(exclusion_log, file.path(data_dir, "study4_exclusions.csv"))

message("Study 4 cleaning complete.")
message("  Staircase long: ", nrow(long),       " rows")
message("  Test phase:     ", nrow(test_out),    " rows")
message("  Summary:        ", nrow(summary_df),  " participants (all retained)")
message("  Flagged:        ", sum(exclusion_log$flag_zero_variance), " zero-variance")
