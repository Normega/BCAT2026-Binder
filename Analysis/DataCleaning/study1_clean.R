# ============================================================
# Study 1 Clean Script — Two-Task Interoception Study
# (co-author's Chapter 3 / formerly Study 2)
#
# Input files (place in base_dir/RawData/):
#   study1_data.xls            — trial-level data with Prolific
#                                demographics already merged (N=184)
#                                (renamed from Study2_Demographics.xls for publication)
#
# Output files (written to data_dir):
#   study1_long.csv        — trial-level rows (Circle1 + Circle2)
#   study1_summary.csv     — one row per participant
#   study1_exclusions.csv  — exclusion log for all 184 participants
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
raw <- read_excel(file.path(raw_dir, "study1_data.xls"),
                  col_types = "text")   # read all as text, coerce below

# Fix leading-space column name on fastORslow
names(raw) <- trimws(names(raw))

# Convert boolean-as-text columns first.
# read_excel() with col_types = "text" reads Excel boolean values as the
# strings "TRUE"/"FALSE"; as.numeric("TRUE") returns NA, so we recode
# these to "1"/"0" before the numeric coercion step below.
bool_cols <- c("DetectACC", "EntrainOK", "BadACC",
               "hardbumper", "easybumper", "ConsentStatus")
for (col in intersect(bool_cols, names(raw))) {
  raw[[col]] <- dplyr::case_when(
    raw[[col]] %in% c("TRUE",  "True",  "true")  ~ "1",
    raw[[col]] %in% c("FALSE", "False", "false") ~ "0",
    TRUE ~ raw[[col]]
  )
}

# Coerce numeric columns
num_cols <- c("Arousal", "Confidence", "DetectACC", "TotalRateChange",
              "TrackACC", "StepSize", "TrialNum", "EventNum",
              "reversalCount", "badcount", "goodcount",
              "EntrainOK", "BadACC",
              paste0("MAIA_", sprintf("%02d", 1:37)),
              "MAIA_total", "age")
for (col in intersect(num_cols, names(raw))) {
  raw[[col]] <- as.numeric(raw[[col]])
}

message("Loaded: ", nrow(raw), " rows, ", n_distinct(raw$UserId), " participants")

# ============================================================
# 2. Exclusion Screening
# ============================================================
# Criterion 1: Zero variance in Arousal ratings across Circle1 trials
#   (non-engagement criterion, consistent across all studies)
circle1_var <- raw |>
  dplyr::filter(Block == "Circle1") |>
  dplyr::group_by(UserId) |>
  dplyr::summarise(
    arousal_sd   = sd(Arousal,    na.rm = TRUE),
    conf_sd      = sd(Confidence, na.rm = TRUE),
    n_c1_trials  = n(),
    .groups = "drop"
  )

# Criterion 2: TrackACC < 70 (no participants fail this in this dataset,
#   but included for transparency and consistency with Study 2)
track_summary <- raw |>
  dplyr::filter(Block == "Circle1") |>
  dplyr::group_by(UserId) |>
  dplyr::summarise(mean_TrackACC = mean(TrackACC, na.rm = TRUE), .groups = "drop")

all_ids <- unique(raw$UserId)

exclusion_log <- tibble(UserId = all_ids) |>
  dplyr::left_join(circle1_var,  by = "UserId") |>
  dplyr::left_join(track_summary, by = "UserId") |>
  dplyr::mutate(
    flag_zero_arousal = !is.na(arousal_sd) & arousal_sd == 0,
    flag_zero_conf    = !is.na(conf_sd)    & conf_sd    == 0,
    flag_low_track    = !is.na(mean_TrackACC) & mean_TrackACC < 70,
    excluded          = flag_zero_arousal | flag_zero_conf | flag_low_track,
    exclusion_reason  = dplyr::case_when(
      flag_zero_arousal ~ "zero_variance_arousal",
      flag_zero_conf    ~ "zero_variance_confidence",
      flag_low_track    ~ "TrackACC_below_70",
      TRUE              ~ "included"
    )
  )

n_excluded <- sum(exclusion_log$excluded)
n_final    <- sum(!exclusion_log$excluded)
message("Excluded: ", n_excluded, " | Final N: ", n_final)
print(exclusion_log |> dplyr::filter(excluded) |> dplyr::select(UserId, exclusion_reason, arousal_sd, mean_TrackACC))

included_ids <- exclusion_log |> dplyr::filter(!excluded) |> dplyr::pull(UserId)

# ============================================================
# 3. Prepare Long-Format Data (Circle1 + Circle2)
# ============================================================
long <- raw |>
  dplyr::filter(Block %in% c("Circle1", "Circle2"),
                UserId %in% included_ids) |>
  dplyr::mutate(
    study_id   = 1L,
    task       = dplyr::recode(Block, "Circle1" = "TaskA", "Circle2" = "TaskB"),
    # Signed change vector: direction × magnitude
    # ChangeType: faster = -1, slower = +1, nochange = 0
    # Note: Study 1 has no salience manipulation — that was introduced in Studies 4 & 5.
    # StepSize here is the adaptive staircase step parameter, not a condition variable.
    NumericChangeType = dplyr::recode(ChangeType,
      "faster"   = -1L,
      "slower"   =  1L,
      "nochange" =  0L
    ),
    ChangeVector = TotalRateChange * NumericChangeType
  ) |>
  dplyr::select(
    study_id, UserId, task, Block, TrialNum, EventNum,
    ChangeType, NumericChangeType, TotalRateChange, ChangeVector,
    DetectACC, Arousal, Confidence, TrackACC, StepSize, reversalCount,
    DetectedChange, DetectedEarly,   # needed for Task B threshold
    age, Sex
  )

message("Long-format rows: ", nrow(long))

# ============================================================
# 4. Compute Summary Variables (one row per participant)
# ============================================================

## 4a. Threshold scores
#
# Task A (Circle1 staircase, 15 trials):
#   Threshold = mean TotalRateChange across the final 4 trials (TrialNum >= 12).
#   This matches prior thesis: "averaging the final 4 values of total breathing
#   rate change in the staircase" (the staircase has exactly 15 trials here).
#
# Task B (Circle2 reaction-time, fixed change magnitude = 0.5):
#   TotalRateChange is always 0.5 and cannot serve as a threshold.
#   DetectedChange records the pulse-rate ratio at the moment the participant
#   pressed the spacebar (final_pulse / initial_pulse). Task B threshold =
#   mean absolute deviation from baseline: mean(|DetectedChange - 1|) across
#   accurately detected change trials (DetectACC == 1, ChangeType != "nochange").
#   Smaller values indicate the participant detected change earlier in its
#   development, i.e., lower threshold / greater sensitivity.

thresh_A <- long |>
  dplyr::filter(task == "TaskA", TrialNum >= 12) |>
  dplyr::group_by(UserId) |>
  dplyr::summarise(thresh_TaskA = mean(TotalRateChange, na.rm = TRUE),
                   .groups = "drop")

thresh_B <- long |>
  dplyr::filter(
    task      == "TaskB",
    DetectACC == 1,
    ChangeType != "nochange"
  ) |>
  dplyr::mutate(DetectedChange = as.numeric(DetectedChange)) |>
  dplyr::group_by(UserId) |>
  dplyr::summarise(
    thresh_TaskB = mean(abs(DetectedChange - 1), na.rm = TRUE),
    .groups = "drop"
  )

thresh <- thresh_A |>
  dplyr::left_join(thresh_B, by = "UserId")

## 4b. Signal Detection Theory (loglinear correction)
sdt <- long |>
  dplyr::filter(ChangeType != "nochange" | !is.na(DetectACC)) |>
  dplyr::mutate(
    signal_trial = ChangeType %in% c("faster", "slower"),
    hit  = signal_trial &  (DetectACC == 1),
    miss = signal_trial & !(DetectACC == 1),
    fa   = !signal_trial & !(DetectACC == 1),
    cr   = !signal_trial &  (DetectACC == 1)
  ) |>
  dplyr::group_by(UserId, task) |>
  dplyr::summarise(
    n_hits   = sum(hit,  na.rm = TRUE),
    n_misses = sum(miss, na.rm = TRUE),
    n_fa     = sum(fa,   na.rm = TRUE),
    n_cr     = sum(cr,   na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    H_corr  = (n_hits   + 0.5) / (n_hits   + n_misses + 1),
    FA_corr = (n_fa     + 0.5) / (n_fa     + n_cr     + 1),
    dprime  = qnorm(H_corr) - qnorm(FA_corr),
    c_bias  = -0.5 * (qnorm(H_corr) + qnorm(FA_corr))
  ) |>
  tidyr::pivot_wider(
    id_cols     = UserId,
    names_from  = task,
    values_from = c(dprime, c_bias, H_corr, FA_corr),
    names_glue  = "{.value}_{task}"
  )

## 4c. Mean Arousal, Confidence, Awareness (cor of Confidence ~ DetectACC per participant-task)
task_summary <- long |>
  dplyr::group_by(UserId, task) |>
  dplyr::summarise(
    mean_Arousal    = mean(Arousal,    na.rm = TRUE),
    mean_Confidence = mean(Confidence, na.rm = TRUE),
    Awareness       = suppressWarnings(cor(Confidence, DetectACC, use = "complete.obs")),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from  = task,
    values_from = c(mean_Arousal, mean_Confidence, Awareness),
    names_glue  = "{.value}_{task}"
  )

## 4d. MAIA subscales (pre-computed in data; pull from questionnaire block)
maia_cols <- c("MAIA_total",
               "MAIA_attentionregulation", "MAIA_bodylisten", "MAIA_emoaware",
               "MAIA_notdistracting", "MAIA_noticing", "MAIA_notworrying",
               "MAIA_selfreg", "MAIA_trusting")
maia <- raw |>
  dplyr::filter(Block == "questionnaire_MAIA", UserId %in% included_ids) |>
  dplyr::group_by(UserId) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(UserId, all_of(intersect(maia_cols, names(raw))))

## 4e. Demographics (one row per participant)
demo <- raw |>
  dplyr::filter(UserId %in% included_ids) |>
  dplyr::group_by(UserId) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(UserId, age, Sex)

## 4f. Join everything into summary
summary_df <- demo |>
  dplyr::left_join(thresh,       by = "UserId") |>
  dplyr::left_join(sdt,          by = "UserId") |>
  dplyr::left_join(task_summary, by = "UserId") |>
  dplyr::left_join(maia,         by = "UserId") |>
  dplyr::mutate(study_id = 1L) |>
  dplyr::select(study_id, UserId, everything())

message("Summary rows: ", nrow(summary_df), " (should be ", n_final, ")")

# ============================================================
# 5. Standardise Output Column Names
# ============================================================
# Renames applied to output files only; internal variable names are unchanged
# so no downstream computation in this script is affected.

long_out <- long |>
  dplyr::rename(
    id       = UserId,
    Accuracy = DetectACC,
    Change   = ChangeVector,
    Group    = task
  ) |>
  dplyr::mutate(
    # Standardise Direction values to match Studies 4 & 5 (title case)
    Direction = dplyr::case_when(
      ChangeType == "faster"   ~ "Faster",
      ChangeType == "slower"   ~ "Slower",
      ChangeType == "nochange" ~ "NoChange",
      TRUE ~ ChangeType
    )
  ) |>
  dplyr::select(-ChangeType)   # Direction replaces ChangeType in output

summary_out <- summary_df |>
  dplyr::rename(
    id                  = UserId,
    # MAIA subscales: lowercase raw names → canonical CamelCase
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
write_csv(long_out,      file.path(data_dir, "study1_long.csv"))
write_csv(summary_out,   file.path(data_dir, "study1_summary.csv"))
write_csv(exclusion_log, file.path(data_dir, "study1_exclusions.csv"))

message("Study 1 cleaning complete.")
message("  Long:       ", nrow(long_out),     " rows")
message("  Summary:    ", nrow(summary_out),  " participants")
message("  Exclusions: ", n_excluded,         " removed (", n_final, " retained)")
