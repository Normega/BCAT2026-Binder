# ============================================================
# Study 5 Clean Script — Fall 2025 In-Person BCAT Study
#
# Input files (place in base_dir/RawData/):
#   study5_longData.csv                  — staircase trial-level data, Sessions 1 & 2
#   study5_testData.csv                  — test-phase trial-level data, Sessions 1 & 2
#   study5_questionnairedata.csv   — reconciled questionnaire composites (N=220)
#
# Note on questionnaire file reconciliation:
#   Two files were provided as study5_questionnairedata:
#     - study5_questionnairedata.xlsx: 0 IDs overlapping with behavioral data;
#       dates 2017-2018. This is the Study 3 (misattribution) questionnaire data
#       and was incorrectly named. Do not use for Study 5.
#     - study5_questionnairedata.csv: 216 of 221 behavioral IDs present. Correct.
#   study5_questionnairedata_clean.csv was produced from the CSV by:
#     - Removing 4 test/dummy entries (IDs 1234, 5678, 99998, 99999)
#     - Standardising column names to match methods doc
#       (Pos→SPANE_Pos, Neg→SPANE_Neg, Stress→BIPS_total, MAIA→MAIA_total, etc.)
#     - Adding has_behavioral_data flag
#   5 quest-only IDs (4604, 13042, 15382, 15724, 16456) are retained in the clean
#   file but will not appear in the final summary (no behavioral data to link to).
#   5 behavioral-only IDs (4606, 11386, 14881, 15559, 16372) have NA questionnaire
#   columns in the summary.
#
# Output files (written to data_dir):
#   study5_long.csv        — staircase trial-level rows (clean sample)
#   study5_test.csv        — test-phase trial-level rows (clean sample)
#   study5_summary.csv     — one row per participant
#   study5_exclusions.csv  — exclusion log for all 221 task participants
#
# Exclusion criteria (per methods doc and pre-registration):
#   1. Pilot participants (research assistants, n=10)
#   2. Test/dummy entries (IDs: 1234, 5678, 99998, 99999; n=1 in task data)
#   3. Zero variance in Arousal ratings across all Session 1 staircase
#      trials (n=3; IDs 6256, 9982, 13711; non-engagement criterion)
#   4. Near-chance accuracy with systematic spamming behaviour (n=1; ID 9553)
#   Final N = 206 behavioural
#
# Notes:
#   - 'VISUAL' recoded to 'visual' in Condition column (capitalisation fix)
#   - Session 2 missingness treated as MAR (session overruns)
#   - Physiological belt data not processed in this script
# ============================================================

# Set Up ---------
packages <- c("readr", "tidyverse", "psych")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
options(readr.show_col_types = FALSE)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

# ============================================================
# 1. Load Raw Data
# ============================================================
staircase  <- read_csv(file.path(raw_dir, "study5_longData.csv"))
test_phase <- read_csv(file.path(raw_dir, "study5_testData.csv"))
quest      <- read_csv(file.path(raw_dir, "study5_questionnairedata.csv"))

# Rename questionnaire columns to canonical cross-study names.
# The rescored CSV uses abbreviated names; map to the standard schema.
quest <- quest |>
  dplyr::rename(
    MAIA_total          = MAIA,
    MAIA_Noticing       = MAIA_Notice,
    MAIA_NotDistracting = MAIA_NotDistract,
    MAIA_NotWorrying    = MAIA_NotWorry,
    MAIA_Trusting       = MAIA_Trust,
    BARQ_total          = BARQ,
    SPANE_Pos           = Pos,
    SPANE_Neg           = Neg,
    BIPS_total          = Stress,
    PHQ4_Anxiety        = Anxiety,
    PHQ4_Depression     = Depression
  )

# Fix capitalisation in Condition
staircase  <- staircase  |> dplyr::mutate(Condition = tolower(Condition))
test_phase <- test_phase |> dplyr::mutate(Condition = tolower(Condition))

message("Staircase:  ", nrow(staircase),  " rows, ", n_distinct(staircase$id),  " participants")
message("Test phase: ", nrow(test_phase),  " rows, ", n_distinct(test_phase$id), " participants")
message("Quest:      ", nrow(quest),       " rows, ", n_distinct(quest$id),      " participants")

# ============================================================
# 2. Exclusion Screening
# ============================================================
all_task_ids <- as.integer(unique(staircase$id))

# Pre-defined exclusion lists (from methods doc and analytical decisions log)
pilot_ids    <- c(5488, 6310, 10105, 10591, 10834, 11890, 12949, 13498, 13879, 15832)
test_ids     <- c(1234, 5678, 99998, 99999)
zero_var_ids <- c(6256, 9982, 13711)
spam_ids     <- c(9553)

# Check for zero-variance in Session 1 staircase Arousal
arousal_var_s1 <- staircase |>
  dplyr::filter(ses == 1) |>
  dplyr::group_by(id) |>
  dplyr::summarise(
    arousal_sd_s1 = sd(as.numeric(Arousal), na.rm = TRUE),
    n_s1_trials   = n(),
    .groups = "drop"
  )

# Compute near-chance accuracy for all participants (ses 1)
acc_s1 <- staircase |>
  dplyr::filter(ses == 1, Direction != "NoChange") |>
  dplyr::group_by(id) |>
  dplyr::summarise(mean_acc_s1 = mean(as.numeric(Accuracy), na.rm = TRUE),
                   .groups = "drop")

# Build exclusion log
exclusion_log <- tibble(id = all_task_ids) |>
  dplyr::left_join(arousal_var_s1, by = "id") |>
  dplyr::left_join(acc_s1,         by = "id") |>
  dplyr::mutate(
    flag_pilot        = id %in% pilot_ids,
    flag_test_entry   = id %in% test_ids,
    flag_zero_arousal = id %in% zero_var_ids,
    flag_near_chance  = id %in% spam_ids,
    # Also catch any additional zero-variance cases not in pre-defined list
    flag_zero_arousal_detected = !is.na(arousal_sd_s1) & arousal_sd_s1 == 0,
    excluded = flag_pilot | flag_test_entry | flag_zero_arousal | flag_near_chance,
    exclusion_reason = dplyr::case_when(
      flag_pilot      ~ "pilot_participant",
      flag_test_entry ~ "test_dummy_entry",
      flag_zero_arousal ~ "zero_variance_arousal",
      flag_near_chance  ~ "near_chance_accuracy_spamming",
      TRUE              ~ "included"
    )
  )

# Flag any zero-variance cases detected algorithmically but not in pre-defined list
additional_zero_var <- exclusion_log |>
  dplyr::filter(flag_zero_arousal_detected & !flag_zero_arousal & !excluded)
if (nrow(additional_zero_var) > 0) {
  warning("Additional zero-variance participants detected but not in pre-defined list: ",
          paste(additional_zero_var$id, collapse = ", "),
          "\nReview these IDs manually.")
}

n_excluded <- sum(exclusion_log$excluded)
n_final    <- sum(!exclusion_log$excluded)
message("\nExclusion summary:")
message("  Pilot participants:       ", sum(exclusion_log$flag_pilot))
message("  Test/dummy entries:       ", sum(exclusion_log$flag_test_entry))
message("  Zero-variance arousal:    ", sum(exclusion_log$flag_zero_arousal))
message("  Near-chance / spamming:   ", sum(exclusion_log$flag_near_chance))
message("  Total excluded:           ", n_excluded)
message("  Final N:                  ", n_final, " (expect 206)")

included_ids <- exclusion_log |> dplyr::filter(!excluded) |> dplyr::pull(id) |> as.integer()

# ============================================================
# 3. Prepare Long-Format Staircase Data
# ============================================================

# Build group lookup first so Condition can be repaired in both long and
# test_out. ses2 rows and ~10 ses1 participants have taskCondition strings
# instead of 'breath'/'visual' in the raw Condition column, which factorise to NA.
group_lookup <- staircase |>
  dplyr::filter(id %in% included_ids, ses == 1) |>
  dplyr::group_by(id) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::transmute(id = as.integer(id),
                   group_ses1 = tolower(Condition))  # "breath" or "visual"

long <- staircase |>
  dplyr::filter(id %in% included_ids) |>
  dplyr::mutate(
    id         = as.integer(id),
    study_id   = 5L,
    level      = as.numeric(level),
    Accuracy   = as.numeric(Accuracy),
    Arousal    = as.numeric(Arousal),
    Confidence = as.numeric(Confidence),
    Change     = as.numeric(Change),
    ses        = factor(ses, levels = c(1, 2), labels = c("ses1", "ses2"))
  ) |>
  dplyr::left_join(group_lookup, by = "id") |>
  dplyr::mutate(
    # Repair Condition: ses2 always breath; ses1 contaminated rows use group_ses1
    Condition = dplyr::case_when(
      tolower(Condition) == "visual"  ~ "visual",
      Condition           == "breath" ~ "breath",
      ses                 == "ses2"   ~ "breath",
      TRUE ~ group_ses1
    ),
    Condition  = factor(Condition, levels = c("breath", "visual")),
    Salience   = factor(Salience,  levels = c("Low", "High")),
    Direction  = factor(Direction, levels = c("Faster", "NoChange", "Slower"))
  ) |>
  dplyr::select(
    study_id, id, ses, taskCondition, Condition, Salience, Direction, DirectionLabel,
    Trial, level, Change, Correct, Accuracy, Response, Arousal, Confidence,
    trial.started, trial.stopped
  )

message("Staircase Condition NA count (expect 0): ",
        sum(is.na(long$Condition)))

# Prepare test-phase data
test_out <- test_phase |>
  dplyr::filter(id %in% included_ids) |>
  dplyr::mutate(id = as.integer(id)) |>
  dplyr::left_join(group_lookup, by = "id") |>
  dplyr::mutate(
    study_id   = 5L,
    level      = as.numeric(level),
    Accuracy   = as.numeric(Accuracy),
    Arousal    = as.numeric(Arousal),
    Confidence = as.numeric(Confidence),
    ses        = factor(ses, levels = c(1, 2), labels = c("ses1", "ses2")),

    # Reconstruct Salience and Direction from taskCondition (fully clean source).
    # The numeric 0/1 and -1/1 columns in the raw file are unreliable; taskCondition
    # is the authoritative encoding for test trials.
    Salience  = factor(
      dplyr::if_else(stringr::str_starts(taskCondition, "high"), "High", "Low"),
      levels = c("Low", "High")
    ),
    Direction = factor(
      dplyr::if_else(stringr::str_ends(taskCondition, "Acc"), "Faster", "Slower"),
      levels = c("Faster", "Slower")
    ),

    # Repair Condition column.
    # Three contamination patterns observed in raw data:
    #   "VISUAL" (uppercase)                          → "visual"
    #   taskCondition strings in ses2 (visual-first)  → "breath" (all ses2 is breath)
    #   taskCondition strings in ses1 (10 participants)→ derive from group_ses1
    Condition = dplyr::case_when(
      tolower(Condition) == "visual"  ~ "visual",
      Condition           == "breath" ~ "breath",
      ses                 == "ses2"   ~ "breath",   # all ses2 is breath
      TRUE ~ group_ses1               # ses1 contaminated rows: use staircase group
    ),
    Condition = factor(Condition, levels = c("breath", "visual"))
  ) |>
  dplyr::select(
    study_id, id, ses, taskCondition, Condition, Salience, Direction, DirectionLabel,
    testTrials.thisN, level, Correct, Accuracy, Response, Arousal, Confidence,
    TestTrial.started, TestTrial.stopped
  )

message("\nStaircase long rows: ", nrow(long))
message("Test phase rows:     ", nrow(test_out))
message("Test Condition NA count (expect 0): ",
        sum(is.na(test_out$Condition)))
message("Test Salience NA count (expect 0):  ",
        sum(is.na(test_out$Salience)))
message("Test Direction NA count (expect 0): ",
        sum(is.na(test_out$Direction)))

# ============================================================
# 4. Compute Summary Variables
# ============================================================

## 4a. Threshold per condition per session (mean level for Trial >= 8)
thresh_staircase <- long |>
  dplyr::filter(as.numeric(as.character(Trial)) >= 8,
                Direction != "NoChange") |>
  dplyr::group_by(id, ses, Salience, Direction) |>
  dplyr::summarise(threshold = mean(level, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(
    names_from  = c(ses, Salience, Direction),
    values_from = threshold,
    names_glue  = "thresh_{ses}_{Salience}_{Direction}"
  )

## 4b. Test-phase threshold per condition per session
thresh_test <- test_out |>
  dplyr::filter(Direction != "NoChange") |>
  dplyr::group_by(id, ses, Salience, Direction) |>
  dplyr::summarise(test_threshold = mean(level, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(
    names_from  = c(ses, Salience, Direction),
    values_from = test_threshold,
    names_glue  = "test_thresh_{ses}_{Salience}_{Direction}"
  )

## 4c. Signal detection per session
sdt <- long |>
  dplyr::mutate(
    signal_trial = Direction != "NoChange",
    hit  = signal_trial  &  (Accuracy == 1),
    miss = signal_trial  & !(Accuracy == 1),
    fa   = !signal_trial & !(Accuracy == 1),
    cr   = !signal_trial &  (Accuracy == 1)
  ) |>
  dplyr::group_by(id, ses) |>
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
  ) |>
  tidyr::pivot_wider(
    id_cols     = id,
    names_from  = ses,
    values_from = c(dprime, c_bias, H_corr, FA_corr),
    names_glue  = "{.value}_{ses}"
  )

## 4d. Mean Arousal, Confidence, Awareness per session
task_summary <- long |>
  dplyr::group_by(id, ses) |>
  dplyr::summarise(
    Condition       = first(as.character(Condition)),
    mean_Arousal    = mean(Arousal,    na.rm = TRUE),
    mean_Confidence = mean(Confidence, na.rm = TRUE),
    Awareness       = suppressWarnings(cor(Confidence, Accuracy, use = "complete.obs")),
    n_trials        = n(),
    n_ses2          = sum(ses == "ses2"),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    id_cols    = id,
    names_from = ses,
    values_from = c(Condition, mean_Arousal, mean_Confidence, Awareness, n_trials),
    names_glue = "{.value}_{ses}"
  )

## 4e. Group assignment (Session 1 condition = group)
group_assign <- long |>
  dplyr::filter(ses == "ses1") |>
  dplyr::group_by(id) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(id, Group = Condition)

## 4f. Session 2 completion flag
ses2_ids <- unique(staircase$id[staircase$ses == 2 & staircase$id %in% included_ids])
group_assign <- group_assign |>
  dplyr::mutate(completed_ses2 = id %in% ses2_ids)

## 4g. Questionnaire data (pre-scored composites; already cleaned)
# Test/dummy entries and quest-only participants (no behavioral data) were
# removed during file reconciliation; included_ids filter handles the rest.
quest_clean <- quest |>
  dplyr::filter(id %in% included_ids) |>
  dplyr::group_by(id) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

## 4h. Join everything
summary_df <- group_assign |>
  dplyr::left_join(task_summary,    by = "id") |>
  dplyr::left_join(thresh_staircase, by = "id") |>
  dplyr::left_join(thresh_test,      by = "id") |>
  dplyr::left_join(sdt,              by = "id") |>
  dplyr::left_join(quest_clean,      by = "id") |>
  dplyr::mutate(study_id = 5L) |>
  dplyr::select(study_id, id, Group, completed_ses2, Age, Gender, everything())

message("Summary rows: ", nrow(summary_df), " (expect 206)")

# ============================================================
# 5. Write Outputs
# ============================================================
write_csv(long,          file.path(data_dir, "study5_long.csv"))
write_csv(test_out,      file.path(data_dir, "study5_test.csv"))
write_csv(summary_df,    file.path(data_dir, "study5_summary.csv"))
write_csv(exclusion_log, file.path(data_dir, "study5_exclusions.csv"))

message("\nStudy 5 cleaning complete.")
message("  Staircase long: ", nrow(long),      " rows")
message("  Test phase:     ", nrow(test_out),   " rows")
message("  Summary:        ", nrow(summary_df), " participants")
message("  Excluded:       ", n_excluded, " removed (", n_final, " retained)")
message("  Session 2 N:    ", sum(summary_df$completed_ses2, na.rm = TRUE))

