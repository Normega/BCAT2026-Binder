# ============================================================
# Study1_MAIA_Reliability.R
# Score the 37-item MAIA and compute Cronbach's alpha
# for total scale and all 8 subscales.
#
# Input:  Study1SourceScales.csv (same directory as script)
# Output: Study1_MAIA_scores.csv     -- one row per participant
#         Study1_MAIA_reliability.csv -- alpha + descriptives per subscale
#
# MAIA SCORING KEY (Mehling et al., 2012; scale 1-5)
# Reverse-scored items (6 - x):
#   Not-Distracting : items 5, 6, 7, 8, 9, 10
#   Not-Worrying    : items 11, 12, 15
# All other items are scored as entered.
#
# Subscale item assignments:
#   Noticing          : 1-4    (4 items)
#   Not-Distracting   : 5-10   (6 items, all reversed)
#   Not-Worrying      : 11-15  (5 items, items 11/12/15 reversed)
#   Attention Reg.    : 16-22  (7 items)
#   Emotional Aware.  : 23-27  (5 items)
#   Self-Regulation   : 28-31  (4 items)
#   Body Listening    : 32-34  (3 items)
#   Trusting          : 35-37  (3 items)
#   Total             : all 37 items (sum)
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
script_dir   <- dirname(rstudioapi::getSourceEditorContext()$path)
input_file   <- file.path(script_dir, "Study1SourceScales.csv")
out_scores   <- file.path(script_dir, "Study1_MAIA_scores.csv")
out_rel      <- file.path(script_dir, "Study1_MAIA_reliability.csv")

# ── Load data ─────────────────────────────────────────────────────────────
raw <- readr::read_csv(input_file)

message(sprintf("Loaded: %d rows, %d participants",
                nrow(raw), dplyr::n_distinct(raw$AuthId)))

# ── Extract one row per participant (MAIA items are constant across trials) ─
# Keep only rows where MAIA_01 is present; one row per AuthId
maia_raw <- raw |>
  dplyr::filter(!is.na(MAIA_01), MAIA_01 != "") |>
  dplyr::distinct(AuthId, .keep_all = TRUE) |>
  dplyr::select(AuthId, dplyr::matches("^MAIA_\\d{2}$")) |>
  dplyr::mutate(dplyr::across(dplyr::matches("^MAIA_\\d{2}$"),
                               as.numeric))

message(sprintf("Participants with MAIA data: %d", nrow(maia_raw)))

# ── Define scoring key ─────────────────────────────────────────────────────
# Items listed as negative strings are reverse-scored (6 - x)
maia_keys <- list(
  MAIA_Noticing     = c("MAIA_01","MAIA_02","MAIA_03","MAIA_04"),
  MAIA_NotDistract  = c("-MAIA_05","-MAIA_06","-MAIA_07",
                        "-MAIA_08","-MAIA_09","-MAIA_10"),
  MAIA_NotWorry     = c("-MAIA_11","-MAIA_12","MAIA_13",
                        "MAIA_14","-MAIA_15"),
  MAIA_AttentionReg = c("MAIA_16","MAIA_17","MAIA_18","MAIA_19",
                        "MAIA_20","MAIA_21","MAIA_22"),
  MAIA_EmoAware     = c("MAIA_23","MAIA_24","MAIA_25","MAIA_26","MAIA_27"),
  MAIA_SelfReg      = c("MAIA_28","MAIA_29","MAIA_30","MAIA_31"),
  MAIA_BodyListen   = c("MAIA_32","MAIA_33","MAIA_34"),
  MAIA_Trusting     = c("MAIA_35","MAIA_36","MAIA_37")
)

# Total key: all items with their reversal status
maia_total_key <- list(
  MAIA_Total = c("-MAIA_05","-MAIA_06","-MAIA_07","-MAIA_08","-MAIA_09","-MAIA_10",
                 "-MAIA_11","-MAIA_12","-MAIA_15",
                 "MAIA_01","MAIA_02","MAIA_03","MAIA_04",
                 "MAIA_13","MAIA_14","MAIA_16","MAIA_17","MAIA_18","MAIA_19",
                 "MAIA_20","MAIA_21","MAIA_22","MAIA_23","MAIA_24","MAIA_25",
                 "MAIA_26","MAIA_27","MAIA_28","MAIA_29","MAIA_30","MAIA_31",
                 "MAIA_32","MAIA_33","MAIA_34","MAIA_35","MAIA_36","MAIA_37")
)

item_data <- maia_raw |>
  dplyr::select(-AuthId) |>
  as.data.frame()

rownames(item_data) <- maia_raw$AuthId

# ── Score subscales (totals, not means) ───────────────────────────────────
# psych::scoreItems with totals = TRUE gives sum scores
subscale_scores <- psych::scoreItems(
  maia_keys, item_data,
  totals  = TRUE,
  min     = 1,
  max     = 5
)

total_scores <- psych::scoreItems(
  maia_total_key, item_data,
  totals = TRUE,
  min    = 1,
  max    = 5
)

scores_df <- tibble::tibble(
  AuthId = maia_raw$AuthId
) |>
  dplyr::bind_cols(
    as.data.frame(subscale_scores$scores),
    as.data.frame(total_scores$scores)
  )

# ── Sanity check against pre-scored values ────────────────────────────────
# Compare MAIA_Total sum to pre-scored MAIA_total column (which is also a sum)
check <- raw |>
  dplyr::filter(!is.na(MAIA_01), MAIA_01 != "") |>
  dplyr::distinct(AuthId, .keep_all = TRUE) |>
  dplyr::select(AuthId, pre_total = MAIA_total) |>
  dplyr::mutate(pre_total = as.numeric(pre_total)) |>
  dplyr::left_join(scores_df |> dplyr::select(AuthId, MAIA_Total), by = "AuthId") |>
  dplyr::mutate(diff = round(abs(pre_total - MAIA_Total), 3))

n_mismatch <- sum(check$diff > 0.01, na.rm = TRUE)
if (n_mismatch == 0) {
  message("Sanity check PASSED: all MAIA_Total values match pre-scored column.")
} else {
  warning(sprintf(
    "Sanity check: %d participants have MAIA_Total discrepancies > 0.01. Check scoring key.",
    n_mismatch
  ))
  print(dplyr::filter(check, diff > 0.01))
}

# ── Cronbach's alpha per subscale ─────────────────────────────────────────
# Apply reversal to item data before passing to psych::alpha
# (psych::alpha handles this via check.keys = TRUE)

all_sub_names <- names(maia_keys)
n_items_map   <- c(
  MAIA_Noticing     = 4L,
  MAIA_NotDistract  = 6L,
  MAIA_NotWorry     = 5L,
  MAIA_AttentionReg = 7L,
  MAIA_EmoAware     = 5L,
  MAIA_SelfReg      = 4L,
  MAIA_BodyListen   = 3L,
  MAIA_Trusting     = 3L,
  MAIA_Total        = 37L
)

.compute_alpha <- function(key_items, item_df, subscale_name) {
  # Strip leading "-" to get item names
  items <- gsub("^-", "", key_items)
  a <- psych::alpha(item_df[, items], check.keys = TRUE)
  vals <- scores_df[[subscale_name]]
  vals <- vals[!is.na(vals)]
  tibble::tibble(
    subscale    = subscale_name,
    n_items     = length(items),
    alpha       = round(a$total$raw_alpha, 3),
    n           = length(vals),
    M           = round(mean(vals), 2),
    SD          = round(sd(vals),   2),
    min_obs     = round(min(vals),  2),
    max_obs     = round(max(vals),  2),
    note        = dplyr::if_else(length(items) <= 3L,
                                  "3-item subscale; alpha interpreted with caution",
                                  "")
  )
}

rel_rows <- dplyr::bind_rows(c(
  lapply(names(maia_keys), function(sub) {
    .compute_alpha(maia_keys[[sub]], item_data, sub)
  }),
  list(.compute_alpha(maia_total_key[["MAIA_Total"]], item_data, "MAIA_Total"))
))

cat("\n=== MAIA Reliability Summary ===\n")
print(as.data.frame(rel_rows |> dplyr::select(subscale, n_items, alpha, n, M, SD)),
      row.names = FALSE)

# ── Save outputs ──────────────────────────────────────────────────────────
readr::write_csv(scores_df, out_scores)
readr::write_csv(rel_rows,  out_rel)

message(sprintf("Saved: %s (%d participants)", basename(out_scores), nrow(scores_df)))
message(sprintf("Saved: %s (%d subscales)",   basename(out_rel),    nrow(rel_rows)))
