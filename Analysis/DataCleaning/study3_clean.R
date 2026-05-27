# ============================================================
# Study 3 Clean Script — Misattribution of Arousal
# (co-author's Chapter 2 / formerly Study 1)
#
# Input files (place in base_dir/RawData/):
#   study3_taskdata.xlsx          — trial-level face rating data (N=103)
#   study3_questionnairedata.xlsx — MAIA-1, PANAS, PSS item responses (N=251)
#
# Output files (written to data_dir):
#   study3_long.csv        — trial-level attraction rating rows
#   study3_summary.csv     — one row per participant
#   study3_exclusions.csv  — exclusion/flag log
#
# Notes:
#   - Participants with only one run (missing "Second Breath Misattribution")
#     are flagged but NOT excluded — they contribute Run 1 data
#   - 4 participants (Pnums 16261, 16283, 17351, 17665) have MA trial data
#     but no questionnaire match; they are retained for trial analyses only
#   - 152 questionnaire-only participants (no MA data) are excluded from
#     all outputs as they have no task data
#   - MAIA scored from 32 individual items (MAIA-1 version, 0-6 scale)
#     Reverse-scored items: 5, 6, 7, 8, 9, 10 (1-indexed)
#
# DATA ISSUE - PANAS item 20 missing:
#   Questionnaires_Clean.xlsx contains only PANAS items 1-19; item 20
#   ("Afraid") is absent from the data file. PANAS subscales are therefore
#   computed from 19 available items only:
#     Positive affect (9 items): 1, 3, 5, 9, 10, 12, 14, 16, 17, 19
#     Negative affect (9 items): 2, 4, 6, 7, 8, 11, 13, 15, 18
#   (Item 20 would normally contribute to the Negative subscale.)
#   Mean-scoring is used so subscale means remain on the 1-5 scale despite
#   the missing item. This omission must be noted in the methods.
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
ma   <- read_excel(file.path(raw_dir, "study3_taskdata.xlsx"))
quest <- read_excel(file.path(raw_dir, "study3_questionnairedata.xlsx"))

# Fix leading-space column name on fastORslow
names(ma) <- trimws(names(ma))

message("MA trial data: ", nrow(ma), " rows, ", n_distinct(ma$Pnum), " participants")
message("Questionnaire data: ", nrow(quest), " rows, ", n_distinct(quest$ID), " participants")

# ============================================================
# 2. Score MAIA-1 (32 items, 0–6 scale) from Questionnaires_Clean
# ============================================================
# MAIA-1 subscale structure (1-indexed item numbers):
#   Noticing:              items  1–4
#   Not-Distracting (R):   items  5–7   (reverse-scored)
#   Not-Worrying (R):      items  8–10  (reverse-scored)
#   Attention Regulation:  items 11–17
#   Emotional Awareness:   items 18–22
#   Self-Regulation:       items 23–25
#   Body Listening:        items 26–28
#   Trusting:              items 29–32
# MAIA-1 items are stored by Qualtrics as 1–7 (representing conceptual 0–6).
# Step 1: recode all 32 items from 1–7 to 0–6 by subtracting 1.
# Step 2: reverse Not-Distracting / Not-Worrying items with 6 - x
#         (correct reversal for a 0–6 scale).
# Without the recode, the old formula 6-x applied to raw 1–7 values produced
# negative scores (e.g., 6 - 7 = -1), corrupting those two subscales.
maia_items    <- paste0("MAIA_", 1:32)
reverse_items <- c(5, 6, 7, 8, 9, 10)  # 1-indexed; Not-Distracting and Not-Worrying

quest_scored <- quest |>
  dplyr::mutate(across(all_of(maia_items), ~ as.numeric(.x) - 1L)) |>   # 1–7 → 0–6
  dplyr::mutate(
    across(all_of(paste0("MAIA_", reverse_items)),
           ~ 6L - .x)                                                    # reverse on 0–6 scale
  ) |>
  dplyr::mutate(
    # Sum scores (0–6 per item after recoding; max per subscale = 6 × n_items)
    MAIA_Noticing       = rowSums(across(paste0("MAIA_", 1:4)),   na.rm = TRUE),  # max 24
    MAIA_NotDistracting = rowSums(across(paste0("MAIA_", 5:7)),   na.rm = TRUE),  # max 18
    MAIA_NotWorrying    = rowSums(across(paste0("MAIA_", 8:10)),  na.rm = TRUE),  # max 18
    MAIA_AttentionReg   = rowSums(across(paste0("MAIA_", 11:17)), na.rm = TRUE),  # max 42
    MAIA_EmoAware       = rowSums(across(paste0("MAIA_", 18:22)), na.rm = TRUE),  # max 30
    MAIA_SelfReg        = rowSums(across(paste0("MAIA_", 23:25)), na.rm = TRUE),  # max 18
    MAIA_BodyListen     = rowSums(across(paste0("MAIA_", 26:28)), na.rm = TRUE),  # max 18
    MAIA_Trusting       = rowSums(across(paste0("MAIA_", 29:32)), na.rm = TRUE),  # max 24
    MAIA_total          = rowSums(across(all_of(maia_items)),     na.rm = TRUE)   # max 192
  )

# Score PANAS from 19 available items (item 20 "Afraid" is absent from the data
# file; see header note above for full documentation of this omission).
# Positive affect (9 items): 1, 3, 5, 9, 10, 12, 14, 16, 17, 19
# Negative affect (9 items): 2, 4, 6, 7, 8, 11, 13, 15, 18   [item 20 omitted]
# Mean-scoring used throughout so values remain on the 1-5 scale.
panas_pos <- c(1, 3, 5, 9, 10, 12, 14, 16, 17, 19)
panas_neg <- c(2, 4, 6, 7,  8, 11, 13, 15, 18)       # item 20 omitted
quest_scored <- quest_scored |>
  dplyr::mutate(across(paste0("PANAS_", c(panas_pos, panas_neg)), as.numeric)) |>
  dplyr::mutate(
    PANAS_Pos   = rowMeans(across(paste0("PANAS_", panas_pos)), na.rm = TRUE),
    PANAS_Neg   = rowMeans(across(paste0("PANAS_", panas_neg)), na.rm = TRUE),
    PANAS_Pos_n = 9L,   # item counts for methods reporting
    PANAS_Neg_n = 9L
  )

# Score PSS-10 (reverse-scored items: 4,5,7,8 using 1-indexed)
pss_reverse <- c(4, 5, 7, 8)
pss_max <- 4L
quest_scored <- quest_scored |>
  dplyr::mutate(across(paste0("PSS", 1:10), as.numeric)) |>
  dplyr::mutate(
    across(paste0("PSS", pss_reverse), ~ pss_max - .x)
  ) |>
  dplyr::mutate(
    PSS_total = rowSums(across(paste0("PSS", 1:10)), na.rm = TRUE)
  )

# ============================================================
# 3. Exclusion / Flag Screening
# ============================================================
all_ma_ids <- unique(ma$Pnum)

# Check run completeness
run_counts <- ma |>
  dplyr::group_by(Pnum) |>
  dplyr::summarise(
    n_runs          = n_distinct(Condition),
    runs_completed  = paste(sort(unique(Condition)), collapse = " | "),
    .groups = "drop"
  )

# Check questionnaire match
quest_ids <- unique(quest$ID)

exclusion_log <- tibble(Pnum = all_ma_ids) |>
  dplyr::left_join(run_counts, by = "Pnum") |>
  dplyr::mutate(
    has_questionnaire  = Pnum %in% quest_ids,
    flag_single_run    = n_runs < 2,
    excluded           = FALSE,   # no exclusions; single-run cases flagged only
    exclusion_reason   = dplyr::case_when(
      flag_single_run & !has_questionnaire ~ "single_run_no_questionnaire",
      flag_single_run                      ~ "single_run_flagged",
      !has_questionnaire                   ~ "no_questionnaire_data",
      TRUE                                 ~ "included"
    )
  )

n_single_run <- sum(exclusion_log$flag_single_run)
n_no_quest   <- sum(!exclusion_log$has_questionnaire)
message("Single-run participants (flagged): ", n_single_run)
message("No questionnaire match: ", n_no_quest)
message("Full-data participants: ", sum(!exclusion_log$flag_single_run & exclusion_log$has_questionnaire))

# ============================================================
# 4. Prepare Long-Format Data
# ============================================================
long <- ma |>
  dplyr::left_join(
    exclusion_log |> dplyr::select(Pnum, flag_single_run, exclusion_reason),
    by = "Pnum"
  ) |>
  dplyr::mutate(
    study_id  = 3L,
    # Signed delta: direction and magnitude of breathing change
    # TotalChange = rate multiplier (< 1 = faster, > 1 = slower, 1 = no change)
    delta     = TotalChange - 1,
    # Direction: derived from sign of delta (TotalChange - 1)
    # delta < 0 = faster breathing, delta > 0 = slower, delta = 0 = no change
    Direction = dplyr::case_when(
      delta < 0  ~ "Faster",
      delta > 0  ~ "Slower",
      delta == 0 ~ "NoChange",
      TRUE       ~ NA_character_
    ),
    # Salience: fastORslow (1 = high/abrupt, 0 = low/gradual)
    Salience  = dplyr::recode(as.character(fastORslow),
                              "1" = "high", "0" = "low"),
    # Run label
    Run = dplyr::recode(Condition,
      "Breath Misattribution"         = "run1",
      "Second Breath Misattribution"  = "run2"
    )
  ) |>
  dplyr::rename(participant_id = Pnum) |>
  dplyr::select(
    study_id, participant_id, Run, Condition, BlockNum,
    gender, SexualOrientation, Salience, Direction, fastORslow,
    TotalChange, delta, FoS_Response, FoS_Accuracy,
    ConfidenceRating, MoodRating, ArousalRating,
    FaceType, FaceRep, Attraction,
    flag_single_run, exclusion_reason
  )

message("Long-format rows: ", nrow(long))

# ============================================================
# 5. Compute Summary Variables
# ============================================================

## 5a. Mean attraction by FaceType and Run
attraction <- long |>
  dplyr::group_by(participant_id, FaceType, Run) |>
  dplyr::summarise(mean_Attraction = mean(Attraction, na.rm = TRUE),
                   .groups = "drop") |>
  tidyr::pivot_wider(
    names_from  = c(FaceType, Run),
    values_from = mean_Attraction,
    names_glue  = "Attraction_{FaceType}_{Run}"
  )

## 5b. Mean Arousal, Mood, and Confidence by run
arousal_summary <- long |>
  dplyr::group_by(participant_id, Run) |>
  dplyr::summarise(
    mean_Arousal     = mean(ArousalRating,    na.rm = TRUE),
    mean_Mood        = mean(MoodRating,        na.rm = TRUE),
    mean_Confidence  = mean(ConfidenceRating,  na.rm = TRUE),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from  = Run,
    values_from = c(mean_Arousal, mean_Mood, mean_Confidence),
    names_glue  = "{.value}_{Run}"
  )

## 5c. Detection accuracy (FoS_Accuracy)
detection <- long |>
  dplyr::group_by(participant_id) |>
  dplyr::summarise(mean_FoS_Accuracy = mean(FoS_Accuracy, na.rm = TRUE),
                   .groups = "drop")

## 5d. Awareness index: within-person cor(ConfidenceRating, FoS_Accuracy)
# Metacognitive accuracy — extent to which trial-level confidence tracks
# detection performance. Computed across all trials (both runs combined).
# Participants with <3 unique values in either variable return NA.
awareness <- long |>
  dplyr::group_by(participant_id) |>
  dplyr::summarise(
    Awareness = suppressWarnings(
      cor(ConfidenceRating, FoS_Accuracy, use = "complete.obs")
    ),
    .groups = "drop"
  )

## 5e. Questionnaire data
quest_out <- quest_scored |>
  dplyr::filter(ID %in% unique(long$participant_id)) |>
  dplyr::group_by(ID) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(ID, Birth_Year, Gender, Orientation, Ethnicity,
                MAIA_total, MAIA_Noticing, MAIA_NotDistracting, MAIA_NotWorrying,
                MAIA_AttentionReg, MAIA_EmoAware, MAIA_SelfReg,
                MAIA_BodyListen, MAIA_Trusting,
                PANAS_Pos, PANAS_Neg, PSS_total, Stress_1) |>
  dplyr::rename(participant_id = ID)

## 5f. Demographics from MA file
demo <- long |>
  dplyr::group_by(participant_id) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(participant_id, gender, SexualOrientation,
                flag_single_run, exclusion_reason)

## 5g. Join
summary_df <- demo |>
  dplyr::left_join(attraction,      by = "participant_id") |>
  dplyr::left_join(arousal_summary, by = "participant_id") |>
  dplyr::left_join(detection,       by = "participant_id") |>
  dplyr::left_join(awareness,       by = "participant_id") |>
  dplyr::left_join(quest_out,       by = "participant_id") |>
  dplyr::mutate(study_id = 3L) |>
  dplyr::select(study_id, participant_id, everything())

message("Summary rows: ", nrow(summary_df))

# ============================================================
# 6. Write Outputs
# ============================================================
# ============================================================
# 6. Standardise Output Column Names
# ============================================================
long_out <- long |>
  dplyr::rename(
    id         = participant_id,
    Change     = delta,
    Accuracy   = FoS_Accuracy,
    Arousal    = ArousalRating,
    Confidence = ConfidenceRating,
    Mood       = MoodRating
  ) |>
  dplyr::mutate(
    Salience = dplyr::case_when(
      Salience == "high" ~ "High",
      Salience == "low"  ~ "Low",
      TRUE ~ Salience
    )
  )

summary_out <- summary_df |>
  dplyr::rename(
    id            = participant_id,
    mean_Accuracy = mean_FoS_Accuracy
  )

# ============================================================
# 7. Write Outputs
# ============================================================
write_csv(long_out,      file.path(data_dir, "study3_long.csv"))
write_csv(summary_out,   file.path(data_dir, "study3_summary.csv"))
write_csv(exclusion_log, file.path(data_dir, "study3_exclusions.csv"))

message("Study 3 cleaning complete.")
message("  Long:         ", nrow(long_out), " rows")
message("  Summary:      ", nrow(summary_df), " participants")
message("  Single-run:   ", n_single_run, " flagged")
message("  No quest:     ", n_no_quest,   " flagged")
