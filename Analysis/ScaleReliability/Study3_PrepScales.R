# ============================================================
# Study3_Scale_Reliability.R
# Score all scales from Study 3 and compute Cronbach's alpha.
#
# Input:  Study3SourceScales.csv (same directory as script)
# Output: Study3_scale_scores.csv    -- one row per participant
#         Study3_scale_reliability.csv -- alpha + descriptives
#
# SCALES:
#   MAIA-1 (32 items, 1-7 scale; Mehling et al., 2012)
#   PSS-10 (10 items, 0-4 scale; Cohen et al., 1983)
#   PANAS  (19 items, 1-5 scale; Watson et al., 1988; item 20 missing)
#   Stress_1 (single VAS item, 0-100; descriptives only)
#
# IMPORTANT NOTES:
#   Study 3 used MAIA-1 (32 items, 1-7 scale). Studies 1-2 used
#   MAIA-2 (37 items, 1-5 scale). Studies 4-5 used Brief MAIA-2
#   (24 items, 0-5 scale). MAIA subscale means should NOT be
#   compared directly across studies without z-score conversion.
#
#   PANAS item 20 ("Afraid") was accidentally omitted in Study 3.
#   The NA subscale therefore has 9 items (not 10). This is noted
#   in reliability output.
#
#   PSS2 column has a trailing space in the source CSV ("PSS2 ").
#   This is handled automatically below.
#
# MAIA-1 SCORING KEY (Mehling et al., 2012; scale 1-7):
#   Reversed items (8-x): 5, 6, 7 (Not-Distracting)
#                          8, 9    (Not-Worrying)
#   Subscales:
#     Noticing          : items 1-4   (4 items)
#     Not-Distracting   : items 5-7   (3 items, all reversed)
#     Not-Worrying      : items 8-10  (5 items, items 8-9 reversed)
#     Attention Reg.    : items 11-17 (7 items)
#     Emotional Aware.  : items 18-22 (5 items)
#     Self-Regulation   : items 23-26 (4 items)
#     Body Listening    : items 27-29 (3 items)
#     Trusting          : items 30-32 (3 items)
#
# PSS-10 SCORING KEY (Cohen et al., 1983; scale 0-4):
#   Reversed items (4-x): 4, 5, 7, 8
#
# PANAS SCORING KEY (Watson et al., 1988; scale 1-5):
#   Positive Affect : items 1, 3, 5, 9, 10, 12, 14, 16, 17, 19
#   Negative Affect : items 2, 4, 6, 7, 8, 11, 13, 15, 18
#                     (item 20 missing; NA subscale has 9 items)
#   No reverse scoring.
# ============================================================

# Set Up ---------
## Load libraries ---------
packages <- c("tidyverse", "psych", "readr")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
options(readr.show_col_types = FALSE)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

# ── Paths ─────────────────────────────────────────────────────────────────
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
input_file <- file.path(script_dir, "Study3SourceScales.csv")
out_scores <- file.path(script_dir, "Study3_scale_scores.csv")
out_rel    <- file.path(script_dir, "Study3_scale_reliability.csv")

# ── Load and clean ─────────────────────────────────────────────────────────
raw <- readr::read_csv(input_file) |>
  dplyr::rename_with(trimws)   # remove trailing whitespace from column names

message(sprintf("Loaded: %d rows, %d participants",
                nrow(raw), dplyr::n_distinct(raw$ID)))

# ── Deduplicate: keep first row per ID ────────────────────────────────────
n_dupes <- sum(duplicated(raw$ID))
if (n_dupes > 0) {
  message(sprintf("Removing %d duplicate ID rows (keeping first occurrence)", n_dupes))
  raw <- raw |> dplyr::distinct(ID, .keep_all = TRUE)
}
message(sprintf("Final sample: %d participants", nrow(raw)))

raw = raw[-which(raw$ID %in% c(11111, 0, 2222, 99999, 77777)),]

# Convert all item columns to numeric
item_cols <- raw |>
  dplyr::select(dplyr::matches("^MAIA_|^PSS|^PANAS_|^Stress_")) |>
  names()

raw <- raw |>
  dplyr::mutate(dplyr::across(dplyr::all_of(item_cols), as.numeric))

# ── MAIA-1 scoring ────────────────────────────────────────────────────────
message("\nScoring MAIA-1...")

maia_items <- paste0("MAIA_", 1:32)

maia_keys <- list(
  MAIA_Noticing     = c("MAIA_1","MAIA_2","MAIA_3","MAIA_4"),
  MAIA_NotDistract  = c("-MAIA_5","-MAIA_6","-MAIA_7"),
  MAIA_NotWorry     = c("-MAIA_8","-MAIA_9","MAIA_10"),
  MAIA_AttentionReg = paste0("MAIA_", 11:17),
  MAIA_EmoAware     = paste0("MAIA_", 18:22),
  MAIA_SelfReg      = paste0("MAIA_", 23:26),
  MAIA_BodyListen   = paste0("MAIA_", 27:29),
  MAIA_Trusting     = paste0("MAIA_", 30:32)
)
maia_total_key <- list(
  MAIA_Total = c(
    "-MAIA_5","-MAIA_6","-MAIA_7",      # Not-Distracting reversed
    "-MAIA_8","-MAIA_9",                  # Not-Worrying reversed
    paste0("MAIA_", c(1:4, 10:32))        # all other items forward
  )
)

maia_item_df <- as.data.frame(dplyr::select(raw, dplyr::all_of(maia_items)))
rownames(maia_item_df) <- raw$ID

maia_sub_scored   <- psych::scoreItems(maia_keys,      maia_item_df,
                                        totals = TRUE, min = 1, max = 7)
maia_total_scored <- psych::scoreItems(maia_total_key, maia_item_df,
                                        totals = TRUE, min = 1, max = 7)

maia_scores <- tibble::tibble(ID = raw$ID) |>
  dplyr::bind_cols(
    as.data.frame(maia_sub_scored$scores),
    as.data.frame(maia_total_scored$scores)
  )

# ── PSS-10 scoring ────────────────────────────────────────────────────────
# Items collected on 1-5 scale, rescored to 0-4 by subtracting 1 per item.
# Reversed items (4, 5, 7, 8) are reversed as 4-x on the 0-4 scale.
# Total range: 0-40.
message("Scoring PSS-10...")

pss_items <- paste0("PSS", 1:10)

# Convert 1-5 to 0-4 before scoring
pss_item_df <- raw |>
  dplyr::select(dplyr::all_of(pss_items)) |>
  dplyr::mutate(dplyr::across(everything(), ~ . - 1)) |>
  as.data.frame()
rownames(pss_item_df) <- raw$ID

pss_keys <- list(
  PSS_total = c("PSS1", "PSS2", "PSS3", "-PSS4", "-PSS5",
                "PSS6", "-PSS7", "-PSS8", "PSS9", "PSS10")
)

pss_scored <- psych::scoreItems(pss_keys, pss_item_df,
                                totals = TRUE, min = 0, max = 4)

pss_scores <- tibble::tibble(
  ID        = raw$ID,
  PSS_total = as.numeric(pss_scored$scores[, "PSS_total"])
)


# ── PANAS scoring ─────────────────────────────────────────────────────────
message("Scoring PANAS...")

panas_pa_items <- paste0("PANAS_", c(1, 3, 5, 9, 10, 12, 14, 16, 17, 19))
panas_na_items <- paste0("PANAS_", c(2, 4, 6, 7, 8, 11, 13, 15, 18))

panas_keys <- list(
  PANAS_PA = panas_pa_items,
  PANAS_NA = panas_na_items
)

panas_all_items <- paste0("PANAS_", 1:19)
panas_item_df   <- as.data.frame(dplyr::select(raw, dplyr::all_of(panas_all_items)))
rownames(panas_item_df) <- raw$ID

panas_scored <- psych::scoreItems(panas_keys, panas_item_df,
                                   totals = TRUE, min = 1, max = 5)

panas_scores <- tibble::tibble(
  ID       = raw$ID,
  PANAS_PA = as.numeric(panas_scored$scores[, "PANAS_PA"]),
  PANAS_NA = as.numeric(panas_scored$scores[, "PANAS_NA"])
)

# ── Stress_1 ──────────────────────────────────────────────────────────────
stress_scores <- tibble::tibble(
  ID       = raw$ID,
  Stress_1 = as.numeric(raw$Stress_1)
)

# ── Combine all scores ────────────────────────────────────────────────────
scores_df <- maia_scores |>
  dplyr::left_join(pss_scores,    by = "ID") |>
  dplyr::left_join(panas_scores,  by = "ID") |>
  dplyr::left_join(stress_scores, by = "ID")

readr::write_csv(scores_df, out_scores)
message(sprintf("Saved: %s (%d participants)", basename(out_scores), nrow(scores_df)))

# ── Alpha computation helper ──────────────────────────────────────────────
.alpha_row <- function(scale, subscale, item_names, item_df,
                        score_col, note = "") {
  sub_df <- item_df[, item_names, drop = FALSE]
  sub_df <- sub_df[complete.cases(sub_df), ]
  if (nrow(sub_df) < 5 || ncol(sub_df) < 2) return(NULL)

  a <- tryCatch(
    psych::alpha(sub_df, check.keys = TRUE),
    error = function(e) { message(sprintf("  [ERROR] %s: %s", subscale, e$message)); NULL }
  )
  if (is.null(a)) return(NULL)

  vals <- score_col[!is.na(score_col)]
  tibble::tibble(
    scale    = scale,
    subscale = subscale,
    n_items  = ncol(sub_df),
    alpha    = round(a$total$raw_alpha, 3),
    n        = length(vals),
    M        = round(mean(vals), 2),
    SD       = round(sd(vals),   2),
    min_obs  = round(min(vals),  2),
    max_obs  = round(max(vals),  2),
    note     = note
  )
}

# ── Compute alpha for all scales ──────────────────────────────────────────
message("\nComputing Cronbach's alpha...")
rel_rows <- list()

# MAIA-1 subscales
for (sub in names(maia_keys)) {
  items <- gsub("^-", "", maia_keys[[sub]])
  rel_rows[[sub]] <- .alpha_row(
    "MAIA-1", sub, items, maia_item_df, scores_df[[sub]],
    note = dplyr::if_else(length(items) <= 3L,
                           "3-item subscale; alpha interpreted with caution", "")
  )
}
# MAIA-1 total
rel_rows[["MAIA_Total"]] <- .alpha_row(
  "MAIA-1", "MAIA_Total",
  gsub("^-","", maia_total_key$MAIA_Total), maia_item_df,
  scores_df$MAIA_Total,
  note = "32-item version; 1-7 scale. Do not compare means with Studies 1-2 (1-5) or 4-5 (0-5)."
)

# PSS-10
pss_item_df_clean <- as.data.frame(dplyr::select(raw, dplyr::all_of(pss_items)))
rel_rows[["PSS"]] <- .alpha_row(
  "PSS-10", "PSS_total", pss_items, pss_item_df_clean, scores_df$PSS_total,
  note = "0-4 scale; items 4,5,7,8 reversed (4-x)"
)

# PANAS PA
rel_rows[["PANAS_PA"]] <- .alpha_row(
  "PANAS", "PANAS_PA", panas_pa_items, panas_item_df, scores_df$PANAS_PA
)
# PANAS NA (9 items -- item 20 missing)
rel_rows[["PANAS_NA"]] <- .alpha_row(
  "PANAS", "PANAS_NA", panas_na_items, panas_item_df, scores_df$PANAS_NA,
  note = "9 items only; item 20 (Afraid) accidentally omitted in Study 3"
)

# Stress_1 -- no alpha (single item)
stress_vals <- scores_df$Stress_1
rel_rows[["Stress_1"]] <- tibble::tibble(
  scale    = "Stress",
  subscale = "Stress_1",
  n_items  = 1L,
  alpha    = NA_real_,
  n        = sum(!is.na(stress_vals)),
  M        = round(mean(stress_vals, na.rm = TRUE), 2),
  SD       = round(sd(stress_vals,   na.rm = TRUE), 2),
  min_obs  = round(min(stress_vals,  na.rm = TRUE), 2),
  max_obs  = round(max(stress_vals,  na.rm = TRUE), 2),
  note     = "Single VAS item (0-100); no alpha"
)

# ── Compile and save ──────────────────────────────────────────────────────
table_reliability <- dplyr::bind_rows(Filter(Negate(is.null), rel_rows))

cat("\n=== Study 3 Scale Reliability Summary ===\n")
print(as.data.frame(
  table_reliability |> dplyr::select(scale, subscale, n_items, alpha, n, M, SD)
), row.names = FALSE)

readr::write_csv(table_reliability, out_rel)
message(sprintf("\nSaved: %s (%d rows)", basename(out_rel), nrow(table_reliability)))
