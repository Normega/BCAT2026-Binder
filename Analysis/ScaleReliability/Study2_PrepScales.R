# ============================================================
# Study2_Scale_Reliability.R
# Parse Study 2 raw JSON, extract all questionnaire scales,
# score MAIA from items, and compute Cronbach's alpha for all.
#
# Input:  Study2raw.json (same directory as script)
# Output: Study2_participants.csv    -- one row per participant,
#                                       all questionnaire items
#                                       + pre-scored values
#         Study2_MAIA_scores.csv     -- MAIA subscale totals
#                                       (re-scored from items)
#         Study2_scale_reliability.csv -- alpha + descriptives
#                                        for all scales
#
# PARSING STRATEGY:
#   Filter to the Block-labelled questionnaire event per
#   participant (e.g. Block == 'questionnaire_MAIA'). Each
#   event carries the full cumulative state, so this gives
#   exactly one clean row per scale per participant.
#
# NOTE on non-MAIA scales:
#   Pre-scored values (REI_EA, BIS, GA7_total, etc.) are
#   used for descriptives as these were computed in the
#   original task script. Cronbach's alpha is computed here
#   from raw items using check.keys = TRUE, which auto-detects
#   reverse-scored items from item-total correlations.
#
# MAIA SCORING KEY (Mehling et al., 2012; scale 1-5):
#   Reversed (6-x): items 5-10 (Not-Distracting);
#                   items 11, 12, 15 (Not-Worrying)
# ============================================================

# Set Up ---------
## Load libraries ---------
packages <- c("tidyverse", "jsonlite", "psych", "readr")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
options(readr.show_col_types = FALSE)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

# ── Paths ─────────────────────────────────────────────────────────────────
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
input_file <- file.path(script_dir, "Study2raw.json")
out_all    <- file.path(script_dir, "Study2_participants.csv")
out_scores <- file.path(script_dir, "Study2_MAIA_scores.csv")
out_rel    <- file.path(script_dir, "Study2_scale_reliability.csv")

# ── Parse JSON ────────────────────────────────────────────────────────────
message("Parsing JSON...")
raw      <- jsonlite::fromJSON(input_file, simplifyVector = FALSE)
sessions <- raw$sessions
message(sprintf("Sessions found: %d", length(sessions)))

# ── Helpers ───────────────────────────────────────────────────────────────

# Find first event matching a Block label
.get_block <- function(events, block_name) {
  for (e in events) {
    if (!is.null(e) && !is.null(e$Block) && e$Block == block_name) return(e)
  }
  NULL
}

# Convert event to one-row tibble for specified fields
.event_to_row <- function(event, fields) {
  vals <- lapply(fields, function(f) {
    if (is.null(event)) return(NA_character_)
    v <- event[[f]]
    if (is.null(v) || identical(v, "")) NA_character_ else as.character(v)
  })
  names(vals) <- fields
  tibble::as_tibble(vals)
}

# ── Field definitions ─────────────────────────────────────────────────────

fields_demo  <- c("Age", "Gender", "Country", "Language",
                  "ConsentId", "ConsentStatus", "StartTime", "UserId")

fields_maia  <- paste0("MAIA_", sprintf("%02d", 1:37))

fields_rei   <- c(
  paste0("EA_",  sprintf("%02d", 1:10)),
  paste0("EE_",  sprintf("%02d", 1:10)),
  paste0("RA_",  sprintf("%02d", 1:10)),
  paste0("RE_",  sprintf("%02d", 1:10)),
  "REI_EA", "REI_EE", "REI_RA", "REI_RE"   # pre-scored means
)

fields_bisbas <- c(
  paste0("BISBAS_", sprintf("%02d", 1:24)),
  "BIS", "BAS_Drive", "BAS_FunSeeking", "BAS_RewardResponsivenes" # pre-scored
)

fields_ga7   <- c(paste0("GA7_0", 1:7), "GA7_total")

fields_phq9  <- c(paste0("PHQ9_0", 1:9), "PHQ9_total")

fields_phq15 <- c(paste0("PHQ15_0", 1:9),
                  paste0("PHQ15_", 10:15), "PHQ15_total")

fields_bart  <- c("BARTtotalpoints", "BARTtotalpumps")

# ── Extract one row per participant ───────────────────────────────────────
message("Extracting questionnaire data...")

participant_rows <- lapply(names(sessions), function(auth_id) {
  events     <- Filter(Negate(is.null), sessions[[auth_id]])

  demo_e   <- .get_block(events, "questionnaire_Demographics")
  maia_e   <- .get_block(events, "questionnaire_MAIA")
  rei_e    <- .get_block(events, "questionnaire_REI")
  bisbas_e <- .get_block(events, "questionnaire_BISBAS")
  ga7_e    <- .get_block(events, "questionnaire_GA7")
  phq9_e   <- .get_block(events, "questionnaire_PHQ9")
  phq15_e  <- .get_block(events, "questionnaire_PHQ15")

  # BART totals from last BART event
  bart_events <- Filter(function(e) identical(e$Block, "BART"), events)
  bart_e <- if (length(bart_events) > 0) bart_events[[length(bart_events)]] else NULL

  tibble::tibble(AuthId = auth_id) |>
    dplyr::bind_cols(
      .event_to_row(demo_e,   fields_demo),
      .event_to_row(maia_e,   fields_maia),
      .event_to_row(rei_e,    fields_rei),
      .event_to_row(bisbas_e, fields_bisbas),
      .event_to_row(ga7_e,    fields_ga7),
      .event_to_row(phq9_e,   fields_phq9),
      .event_to_row(phq15_e,  fields_phq15),
      .event_to_row(bart_e,   fields_bart)
    )
})

participants <- dplyr::bind_rows(participant_rows)

# Convert numeric columns
participants <- participants |>
  dplyr::mutate(dplyr::across(
    -c(AuthId, Gender, Country, Language, ConsentId, ConsentStatus, UserId),
    ~ suppressWarnings(as.numeric(.))
  ))

n_maia <- sum(!is.na(participants$MAIA_01))
message(sprintf("Participants: %d total, %d with MAIA data", nrow(participants), n_maia))

readr::write_csv(participants, out_all)
message(sprintf("Saved: %s", basename(out_all)))

# ── MAIA: re-score from items ─────────────────────────────────────────────
maia_df <- participants |>
  dplyr::filter(!is.na(MAIA_01)) |>
  dplyr::select(AuthId, dplyr::all_of(fields_maia))

maia_items <- as.data.frame(dplyr::select(maia_df, -AuthId))
rownames(maia_items) <- maia_df$AuthId

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
maia_total_key <- list(
  MAIA_Total = c(
    "-MAIA_05","-MAIA_06","-MAIA_07","-MAIA_08","-MAIA_09","-MAIA_10",
    "-MAIA_11","-MAIA_12","-MAIA_15",
    "MAIA_01","MAIA_02","MAIA_03","MAIA_04",
    "MAIA_13","MAIA_14","MAIA_16","MAIA_17","MAIA_18","MAIA_19",
    "MAIA_20","MAIA_21","MAIA_22","MAIA_23","MAIA_24","MAIA_25",
    "MAIA_26","MAIA_27","MAIA_28","MAIA_29","MAIA_30","MAIA_31",
    "MAIA_32","MAIA_33","MAIA_34","MAIA_35","MAIA_36","MAIA_37"
  )
)

maia_sub_scored   <- psych::scoreItems(maia_keys,     maia_items, totals = TRUE, min = 1, max = 5)
maia_total_scored <- psych::scoreItems(maia_total_key, maia_items, totals = TRUE, min = 1, max = 5)

scores_df <- tibble::tibble(AuthId = maia_df$AuthId) |>
  dplyr::bind_cols(
    as.data.frame(maia_sub_scored$scores),
    as.data.frame(maia_total_scored$scores)
  )

readr::write_csv(scores_df, out_scores)
message(sprintf("Saved: %s (%d participants)", basename(out_scores), nrow(scores_df)))

# ── Alpha computation helper ──────────────────────────────────────────────
# Computes Cronbach's alpha from a matrix of raw item responses.
# Uses check.keys = TRUE to auto-detect and flip negatively keyed items.
# For descriptives, uses the pre-scored column (pre_col) from participants.
.alpha_row <- function(scale, subscale, item_cols, item_df,
                       pre_col = NULL, n_items = NULL, note = "") {
  # Keep only rows where all items are non-NA
  sub_df <- item_df[, item_cols, drop = FALSE]
  sub_df <- sub_df[complete.cases(sub_df), ]
  if (nrow(sub_df) < 5 || ncol(sub_df) < 2) {
    message(sprintf("  [SKIP] %s %s: too few complete cases or items", scale, subscale))
    return(NULL)
  }

  a <- tryCatch(
    psych::alpha(sub_df, check.keys = TRUE),
    error = function(e) { message(sprintf("  [ERROR] %s %s: %s", scale, subscale, e$message)); NULL }
  )
  if (is.null(a)) return(NULL)

  # Descriptives from pre-scored column if provided, else from scored item means
  if (!is.null(pre_col) && !all(is.na(pre_col))) {
    vals <- pre_col[!is.na(pre_col)]
  } else {
    vals <- rowSums(sub_df)
  }

  n_items_val <- if (!is.null(n_items)) n_items else ncol(sub_df)

  tibble::tibble(
    scale      = scale,
    subscale   = subscale,
    n_items    = n_items_val,
    alpha      = round(a$total$raw_alpha, 3),
    n          = length(vals),
    M          = round(mean(vals), 2),
    SD         = round(sd(vals),   2),
    min_obs    = round(min(vals),  2),
    max_obs    = round(max(vals),  2),
    note       = note
  )
}

# ── Compute alpha for all scales ──────────────────────────────────────────
message("\nComputing Cronbach's alpha for all scales...")

all_items <- as.data.frame(
  participants |> dplyr::select(-AuthId, -Gender, -Country,
                                 -Language, -ConsentId, -ConsentStatus, -UserId)
)

rel_rows <- list()

# ── MAIA (known key, totals) ──────────────────────────────────────────────
for (sub in names(maia_keys)) {
  items <- gsub("^-", "", maia_keys[[sub]])
  a     <- psych::alpha(maia_items[, items], check.keys = TRUE)
  vals  <- scores_df[[sub]]
  vals  <- vals[!is.na(vals)]
  rel_rows[[sub]] <- tibble::tibble(
    scale    = "MAIA",
    subscale = sub,
    n_items  = length(items),
    alpha    = round(a$total$raw_alpha, 3),
    n        = length(vals),
    M        = round(mean(vals), 2),
    SD       = round(sd(vals),   2),
    min_obs  = round(min(vals),  2),
    max_obs  = round(max(vals),  2),
    note     = dplyr::if_else(length(items) <= 3L,
                               "3-item subscale; alpha interpreted with caution", "")
  )
}
# MAIA total
a_total <- psych::alpha(maia_items[, gsub("^-","",maia_total_key$MAIA_Total)], check.keys = TRUE)
vals_t  <- scores_df$MAIA_Total; vals_t <- vals_t[!is.na(vals_t)]
rel_rows[["MAIA_Total"]] <- tibble::tibble(
  scale="MAIA", subscale="MAIA_Total", n_items=37L,
  alpha=round(a_total$total$raw_alpha,3), n=length(vals_t),
  M=round(mean(vals_t),2), SD=round(sd(vals_t),2),
  min_obs=round(min(vals_t),2), max_obs=round(max(vals_t),2), note=""
)

# ── REI subscales ─────────────────────────────────────────────────────────
# Items are scored 1-5; inv_ columns are pre-reversed versions (not re-used here)
rei_subs <- list(
  REI_EA = paste0("EA_", sprintf("%02d", 1:10)),
  REI_EE = paste0("EE_", sprintf("%02d", 1:10)),
  REI_RA = paste0("RA_", sprintf("%02d", 1:10)),
  REI_RE = paste0("RE_", sprintf("%02d", 1:10))
)
for (sub in names(rei_subs)) {
  rel_rows[[sub]] <- .alpha_row(
    "REI", sub, rei_subs[[sub]], all_items,
    pre_col = participants[[sub]],
    note = "alpha via check.keys; pre-scored means used for descriptives"
  )
}

# ── BIS/BAS subscales ─────────────────────────────────────────────────────
# Carver & White (1994): fillers = items 1, 6, 11, 17 (excluded)
# BIS: 2, 8, 13, 16, 19, 22, 24 (items 2, 22 reverse-scored)
# BAS Drive: 3, 9, 12, 21
# BAS Fun Seeking: 5, 10, 15, 20
# BAS Reward Resp: 4, 7, 14, 18, 23
bisbas_subs <- list(
  BIS = paste0("BISBAS_", sprintf("%02d", c(2, 8, 13, 16, 19, 22, 24))),
  BAS_Drive    = paste0("BISBAS_", sprintf("%02d", c(3, 9, 12, 21))),
  BAS_FunSeeking = paste0("BISBAS_", sprintf("%02d", c(5, 10, 15, 20))),
  BAS_RewardResp = paste0("BISBAS_", sprintf("%02d", c(4, 7, 14, 18, 23)))
)
bisbas_pre <- list(
  BIS            = "BIS",
  BAS_Drive      = "BAS_Drive",
  BAS_FunSeeking = "BAS_FunSeeking",
  BAS_RewardResp = "BAS_RewardResponsivenes"
)
for (sub in names(bisbas_subs)) {
  rel_rows[[paste0("BISBAS_", sub)]] <- .alpha_row(
    "BIS/BAS", sub, bisbas_subs[[sub]], all_items,
    pre_col = participants[[bisbas_pre[[sub]]]],
    note = "alpha via check.keys; pre-scored means used for descriptives"
  )
}

# ── GAD-7 ─────────────────────────────────────────────────────────────────
rel_rows[["GAD7"]] <- .alpha_row(
  "GAD-7", "GAD7_total",
  paste0("GA7_0", 1:7), all_items,
  pre_col = participants$GA7_total,
  n_items = 7L
)

# ── PHQ-9 ─────────────────────────────────────────────────────────────────
rel_rows[["PHQ9"]] <- .alpha_row(
  "PHQ-9", "PHQ9_total",
  paste0("PHQ9_0", 1:9), all_items,
  pre_col = participants$PHQ9_total,
  n_items = 9L
)

# ── PHQ-15 ────────────────────────────────────────────────────────────────
# PHQ15_total from JSON is corrupted for 2 participants (missing item 04
# produced garbage pre-scored values). Re-score from raw items instead.
# Participants with any missing item receive NA (no partial scoring).
phq15_item_df <- participants |>
  dplyr::select(dplyr::all_of(phq15_items)) |>
  dplyr::mutate(dplyr::across(everything(), as.numeric))

phq15_total_computed <- rowSums(phq15_item_df, na.rm = FALSE)

cat(sprintf(
  "PHQ-15 re-scored from items: %d complete, %d with missing items (set to NA)\n",
  sum(!is.na(phq15_total_computed)),
  sum(is.na(phq15_total_computed))
))

phq15_items <- c(paste0("PHQ15_0", 1:9), paste0("PHQ15_", 10:15))
rel_rows[["PHQ15"]] <- .alpha_row(
  "PHQ-15", "PHQ15_total",
  phq15_items, all_items,
  pre_col = phq15_total_computed,
  n_items = 15L,
  note = "Re-scored from items; 2 participants excluded due to missing PHQ15_04"
)

# ── BART (no alpha -- behavioral task) ───────────────────────────────────
bart_vals_pts  <- participants$BARTtotalpoints
bart_vals_pmp  <- participants$BARTtotalpumps
rel_rows[["BART_points"]] <- tibble::tibble(
  scale="BART", subscale="BARTtotalpoints", n_items=NA_integer_,
  alpha=NA_real_, n=sum(!is.na(bart_vals_pts)),
  M=round(mean(bart_vals_pts,na.rm=TRUE),2), SD=round(sd(bart_vals_pts,na.rm=TRUE),2),
  min_obs=round(min(bart_vals_pts,na.rm=TRUE),2), max_obs=round(max(bart_vals_pts,na.rm=TRUE),2),
  note="behavioral task; no alpha"
)
rel_rows[["BART_pumps"]] <- tibble::tibble(
  scale="BART", subscale="BARTtotalpumps", n_items=NA_integer_,
  alpha=NA_real_, n=sum(!is.na(bart_vals_pmp)),
  M=round(mean(bart_vals_pmp,na.rm=TRUE),2), SD=round(sd(bart_vals_pmp,na.rm=TRUE),2),
  min_obs=round(min(bart_vals_pmp,na.rm=TRUE),2), max_obs=round(max(bart_vals_pmp,na.rm=TRUE),2),
  note="behavioral task; no alpha"
)

# ── Compile and save ──────────────────────────────────────────────────────
table_reliability <- dplyr::bind_rows(Filter(Negate(is.null), rel_rows))

cat("\n=== Study 2 Scale Reliability Summary ===\n")
print(as.data.frame(
  table_reliability |> dplyr::select(scale, subscale, n_items, alpha, n, M, SD)
), row.names = FALSE)

readr::write_csv(table_reliability, out_rel)
message(sprintf("\nSaved: %s (%d rows)", basename(out_rel), nrow(table_reliability)))
