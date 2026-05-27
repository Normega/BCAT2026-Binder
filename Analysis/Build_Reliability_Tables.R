# ============================================================
# Build_Reliability_Tables.R  (v2 - APA7 style)
# Reads all scale reliability CSVs and generates two reliability
# tables as supplementary Word document.
#
# Changes from v1:
#   - APA7 three-rule style (white, no colour blocks)
#   - Removed from Table SR2: BIS/BAS, BART, GAD-7, PHQ-9,
#     PHQ-15, REI (per SM2: "administered as exploratory
#     individual difference measures and are not included
#     in the main reliability summary")
#
# Input:  Scale reliability CSVs in ScaleReliability directory
# Output: Scale_Reliability_Tables.docx
# ============================================================

# Set Up ---------

# ── Config ────────────────────────────────────────────────────────────────────
FONT       <- "Times New Roman"
SZ         <- 11    # slightly smaller for dense reliability tables
SZ_NOTE    <- 9

script_dir <- file.path(BASE_DIR, "Data", "ScaleReliability")
output_dir  <- file.path(BASE_DIR, "Tables")

NS <- c(s1 = 184, s2 = 167, s3 = 246, s4 = 146, s5 = 225)

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

# ── APA7 style ────────────────────────────────────────────────────────────────
apa7_style <- function(ft) {
  ft |>
    font(fontname = FONT, part = "all") |>
    fontsize(size = SZ, part = "all") |>
    bg(bg = "white", part = "all") |>
    color(color = "black", part = "all") |>
    border_remove() |>
    hline_top(border = officer::fp_border(width = 1.5, color = "black"), part = "header") |>
    hline_bottom(border = officer::fp_border(width = 0.75, color = "black"), part = "header") |>
    hline_bottom(border = officer::fp_border(width = 1.5, color = "black"), part = "body") |>
    align(align = "left",   j = 1, part = "all") |>
    align(align = "center", j = -1, part = "all") |>
    padding(padding.top = 3, padding.bottom = 3,
            padding.left = 5, padding.right = 5, part = "all")
}

# ── Helpers ───────────────────────────────────────────────────────────────────
fmt_alpha <- function(x) {
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v) || v <= 0) return("--")
  sub("^0", "", sprintf("%.2f", v))
}

fmt_r <- function(x) {
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v) || v <= 0) return("--")
  paste0("r\u00a0=\u00a0", sub("^0", "", sprintf("%.2f", v)), "\u2020")
}

load_rel <- function(filename) {
  fp <- file.path(script_dir, filename)
  if (!file.exists(fp)) {
    warning(sprintf("File not found: %s", filename))
    return(tibble::tibble())
  }
  df <- readr::read_csv(fp) |>
    dplyr::mutate(dplyr::across(everything(), as.character))
  if (!"scale"       %in% names(df)) df$scale       <- NA_character_
  if (!"interitem_r" %in% names(df)) df$interitem_r <- NA_character_
  df
}

s1 <- load_rel("Study1_scale_reliability.csv")
s2 <- load_rel("Study2_scale_reliability.csv")
s3 <- load_rel("Study3_scale_reliability.csv")
s4 <- load_rel("Study4_scale_reliability.csv")
s5 <- load_rel("Study5_scale_reliability.csv")

get_val <- function(study_df, sub_key) {
  if (nrow(study_df) == 0) return("--")
  sub_key <- sub_key[!is.na(sub_key) & sub_key != ""]
  if (length(sub_key) == 0) return("--")

  row <- study_df |>
    dplyr::filter(subscale %in% sub_key) |>
    dplyr::slice(1)
  if (nrow(row) == 0) return("--")

  note_val <- row$note %||% ""
  if (!is.na(note_val) && grepl("behavioral|Single VAS", note_val, ignore.case = TRUE)) {
    return("n/a")
  }
  ni <- suppressWarnings(as.integer(row$n_items))
  ir <- row$interitem_r %||% NA_character_
  if ((!is.na(ni) && ni == 2L) ||
      (!is.na(ir) && ir != "" && ir != "NA")) {
    return(fmt_r(ir))
  }
  fmt_alpha(row$alpha)
}

# ══════════════════════════════════════════════════════════════════════════════
# TABLE SR1: MAIA Internal Consistency
# ══════════════════════════════════════════════════════════════════════════════

maia_subs <- list(
  list(d = "Noticing",             n2 = 4,  n1 = 4,  nb = 3,
       s1 = "MAIA_Noticing",    s2 = "MAIA_Noticing",    s3 = "MAIA_Noticing",
       s4 = "MAIA_Noticing",    s5 = "MAIA_Notice"),
  list(d = "Not-Distracting",      n2 = 6,  n1 = 3,  nb = 3,
       s1 = "MAIA_NotDistract", s2 = "MAIA_NotDistract", s3 = "MAIA_NotDistract",
       s4 = "MAIA_NotDistract", s5 = "MAIA_NotDistract"),
  list(d = "Not-Worrying",         n2 = 5,  n1 = 3,  nb = 3,
       s1 = "MAIA_NotWorry",    s2 = "MAIA_NotWorry",    s3 = "MAIA_NotWorry",
       s4 = "MAIA_NotWorry",    s5 = "MAIA_NotWorry"),
  list(d = "Attention Regulation", n2 = 7,  n1 = 7,  nb = 3,
       s1 = "MAIA_AttentionReg", s2 = "MAIA_AttentionReg", s3 = "MAIA_AttentionReg",
       s4 = "MAIA_AttentionReg", s5 = "MAIA_AttentionReg"),
  list(d = "Emotional Awareness",  n2 = 5,  n1 = 5,  nb = 3,
       s1 = "MAIA_EmoAware",    s2 = "MAIA_EmoAware",    s3 = "MAIA_EmoAware",
       s4 = "MAIA_EmoAware",    s5 = "MAIA_EmoAware"),
  list(d = "Self-Regulation",      n2 = 4,  n1 = 4,  nb = 3,
       s1 = "MAIA_SelfReg",     s2 = "MAIA_SelfReg",     s3 = "MAIA_SelfReg",
       s4 = "MAIA_SelfReg",     s5 = "MAIA_SelfReg"),
  list(d = "Body Listening",       n2 = 3,  n1 = 3,  nb = 3,
       s1 = "MAIA_BodyListen",  s2 = "MAIA_BodyListen",  s3 = "MAIA_BodyListen",
       s4 = "MAIA_BodyListen",  s5 = "MAIA_BodyListen"),
  list(d = "Trusting",             n2 = 3,  n1 = 3,  nb = 3,
       s1 = "MAIA_Trusting",    s2 = "MAIA_Trusting",    s3 = "MAIA_Trusting",
       s4 = "MAIA_Trust",       s5 = "MAIA_Trust"),
  list(d = "Total",                n2 = 37, n1 = 32, nb = 24,
       s1 = "MAIA_Total",       s2 = "MAIA_Total",       s3 = "MAIA_Total",
       s4 = "MAIA",             s5 = "MAIA",
       is_total = TRUE)
)

build_sr1 <- function() {
  rows        <- list()
  header_rows <- integer(0)
  total_rows  <- integer(0)
  ri <- 0L

  add_ver_header <- function(label) {
    ri <<- ri + 1L
    header_rows <<- c(header_rows, ri)
    rows[[ri]] <<- tibble::tibble(
      subscale = label, items = "",
      s1 = "", s2 = "", s3 = "", s4 = "", s5 = ""
    )
  }

  add_row <- function(sub, items_val, v1, v2, v3, v4, v5) {
    ri <<- ri + 1L
    if (isTRUE(sub$is_total)) total_rows <<- c(total_rows, ri)
    rows[[ri]] <<- tibble::tibble(
      subscale = sub$d, items = as.character(items_val),
      s1 = v1, s2 = v2, s3 = v3, s4 = v4, s5 = v5
    )
  }

  # MAIA-2 (37 items; Studies 1-2)
  add_ver_header("MAIA-2 (37 items)   Mehling et al., 2018; scale 0\u20135; Studies 1\u20132")
  for (sub in maia_subs) {
    add_row(sub, sub$n2,
            get_val(s1, sub$s1), get_val(s2, sub$s2), "--", "--", "--")
  }

  # MAIA-1 (32 items; Study 3)
  add_ver_header("MAIA-1 (32 items)   Mehling et al., 2012; scale 1\u20137 in Study 3; means not comparable with MAIA-2")
  for (sub in maia_subs) {
    add_row(sub, sub$n1, "--", "--", get_val(s3, sub$s3), "--", "--")
  }

  # Brief MAIA-2 (24 items; Studies 4-5)
  add_ver_header("Brief MAIA-2 (24 items)   Rogowska et al., 2023; 3 items per subscale; scale 0\u20135; Studies 4\u20135")
  for (sub in maia_subs) {
    add_row(sub, sub$nb, "--", "--", "--",
            get_val(s4, sub$s4 %||% NA),
            get_val(s5, sub$s5))
  }

  list(data = dplyr::bind_rows(rows),
       header_rows = header_rows,
       total_rows  = total_rows)
}

make_sr1 <- function() {
  obj <- build_sr1()
  df  <- obj$data
  hr  <- obj$header_rows
  tr  <- obj$total_rows
  n   <- nrow(df)

  ft <- flextable(df) |>
    set_header_labels(
      subscale = "Subscale", items = "Items",
      s1 = sprintf("Study 1\n(N\u202f=\u202f%d)", NS["s1"]),
      s2 = sprintf("Study 2\n(N\u202f=\u202f%d)", NS["s2"]),
      s3 = sprintf("Study 3\n(N\u202f=\u202f%d)", NS["s3"]),
      s4 = sprintf("Study 4\n(N\u202f=\u202f%d)", NS["s4"]),
      s5 = sprintf("Study 5\n(N\u202f=\u202f%d)", NS["s5"])
    ) |>
    # Version header rows: merge across all columns, bold, indent
    merge_h(i = hr) |>
    bold(i = hr) |>
    italic(i = hr) |>
    # Total rows: bold
    bold(i = tr) |>
    # Column widths (landscape 9.5" content)
    width(j = "subscale", width = 2.5) |>
    width(j = "items",    width = 0.4) |>
    width(j = c("s1","s2","s3","s4","s5"), width = 1.32) |>
    apa7_style()

  # Add thin horizontal rules before each version header (except first)
  for (h in hr[-1]) {
    ft <- hline(ft, i = h - 1,
                border = officer::fp_border(width = 0.75, color = "black"),
                part = "body")
  }
  ft
}

# ══════════════════════════════════════════════════════════════════════════════
# TABLE SR2: All Other Measures
# ── Omitted per SM2: BIS/BAS, BART, GAD-7, PHQ-9, PHQ-15, REI ──────────────
# ══════════════════════════════════════════════════════════════════════════════

other_rows_def <- list(
  # BIPS
  list("BIPS",     "Total",              "9",          NA,           NA,          "Stress",        "Stress"),
  # BARQ-R
  list("BARQ-R",   "Total",              "12",         NA,           NA,          "BARQ",          "BARQ"),
  # MSES-12
  list("MSES-12",  "Total",              "12",         NA,           NA,          NA,              "MSES"),
  list("",         "Self-Doubt",         "2",          NA,           NA,          NA,              "MSES_selfdoubt"),
  # PANAS
  list("PANAS",    "Positive Affect",    "10",         NA,           "PANAS_PA",  NA,              NA),
  list("",         "Negative Affect",    "9\u2021",    NA,           "PANAS_NA",  NA,              NA),
  # PAQ-S
  list("PAQ-S",    "Total",              "6",          NA,           NA,          NA,              "Alexithymia"),
  # PHQ-4
  list("PHQ-4",    "Anxiety",            "2",          NA,           NA,          "Anxiety",       "PHQ4_Anxiety"),
  list("",         "Depression",         "2",          NA,           NA,          "Depression",    "PHQ4_Depression"),
  # PSS-10
  list("PSS-10",   "Total",              "10",         NA,           "PSS_total", NA,              NA),
  # SPANE
  list("SPANE",    "Positive",           "6",          NA,           NA,          "Pos",           "Pos"),
  list("",         "Negative",           "6",          NA,           NA,          "Neg",           "Neg"),
  # Stress (Study 3 single VAS item)
  list("Stress",   "Single item (VAS)",  "1",          NA,           "Stress_1",  NA,              NA),
  # Wellbeing (Study 4 composite)
  list("Wellbeing","Total",              "8",          NA,           NA,          "Wellbeing",     NA)
)

build_sr2 <- function() {
  purrr::map_dfr(other_rows_def, function(r) {
    tibble::tibble(
      scale    = as.character(r[[1]]),
      subscale = as.character(r[[2]]),
      items    = as.character(r[[3]]),
      s2       = get_val(s2, r[[4]]),
      s3       = get_val(s3, r[[5]]),
      s4       = get_val(s4, r[[6]]),
      s5       = get_val(s5, r[[7]])
    )
  })
}

make_sr2 <- function() {
  df <- build_sr2()

  flextable(df) |>
    set_header_labels(
      scale = "Scale", subscale = "Subscale", items = "Items",
      s2 = sprintf("Study 2\n(N\u202f=\u202f%d)", NS["s2"]),
      s3 = sprintf("Study 3\n(N\u202f=\u202f%d)", NS["s3"]),
      s4 = sprintf("Study 4\n(N\u202f=\u202f%d)", NS["s4"]),
      s5 = sprintf("Study 5\n(N\u202f=\u202f%d)", NS["s5"])
    ) |>
    merge_v(j = "scale") |>
    valign(j = "scale", valign = "top", part = "body") |>
    bold(i = which(df$scale != ""), j = "scale") |>
    width(j = "scale",    width = 0.9) |>
    width(j = "subscale", width = 2.0) |>
    width(j = "items",    width = 0.45) |>
    width(j = c("s2","s3","s4","s5"), width = 1.54) |>
    apa7_style()
}

# ── Build Word document ───────────────────────────────────────────────────────
message("Building reliability tables Word document...")

ps_landscape <- officer::prop_section(
  page_size    = officer::page_size(width = 11, height = 8.5, orient = "landscape"),
  page_margins = officer::page_mar(top = 0.75, bottom = 0.75, left = 0.75, right = 0.75)
)

note_prop <- officer::fp_text(font.family = FONT, font.size = SZ_NOTE)
bold_prop <- officer::fp_text(font.family = FONT, font.size = SZ_NOTE, bold = TRUE)

doc <- officer::read_docx() |>
  officer::body_set_default_section(ps_landscape) |>

  officer::body_add_par("Table SR1", style = "heading 3") |>
  officer::body_add_par(
    "MAIA Internal Consistency Across Studies and Versions",
    style = "Normal"
  ) |>
  flextable::body_add_flextable(make_sr1()) |>
  officer::body_add_par("", style = "Normal") |>
  officer::body_add_fpar(officer::fpar(
    officer::ftext("Note. ", prop = bold_prop),
    officer::ftext(paste0(
      "\u03b1 = Cronbach\u2019s alpha.  -- = scale not administered.  ",
      "MAIA-1: original 32-item version (Mehling et al., 2012).  ",
      "MAIA-2: revised 37-item version (Mehling et al., 2018).  ",
      "Brief MAIA-2: 24-item version with 3 items per subscale (Rogowska et al., 2023).  ",
      "Subscale means are not comparable across versions (MAIA-1: 1\u20137; MAIA-2: 0\u20135).  ",
      "3-item subscales (Brief MAIA-2) should be interpreted with caution."
    ), prop = note_prop)
  )) |>

  officer::body_add_par("", style = "Normal") |>

  officer::body_add_par("Table SR2", style = "heading 3") |>
  officer::body_add_par(
    "Scale Internal Consistency: All Other Measures",
    style = "Normal"
  ) |>
  flextable::body_add_flextable(make_sr2()) |>
  officer::body_add_par("", style = "Normal") |>
  officer::body_add_fpar(officer::fpar(
    officer::ftext("Note. ", prop = bold_prop),
    officer::ftext(paste0(
      "\u03b1 = Cronbach\u2019s alpha.  -- = scale not administered.  ",
      "n/a = not applicable (behavioral task or single item).  ",
      "\u2020 Inter-item Pearson r reported; \u03b1 not appropriate for 2-item scales.  ",
      "\u2021 PANAS Negative Affect has 9 items in Study 3; item 20 (\u201cAfraid\u201d) was inadvertently omitted.  ",
      "BIPS = Brief Inventory of Perceived Stress; BARQ-R = Body Awareness Rating Questionnaire\u2013Revised; ",
      "PAQ-S = Perth Alexithymia Questionnaire\u2013Short Form; MSES-12 = Multidimensional Self-Esteem Scale ",
      "ultra-short form; SPANE = Scale of Positive and Negative Experience.  ",
      "BIS/BAS, BART, GAD-7, PHQ-9, PHQ-15, and REI were administered as exploratory individual ",
      "difference measures and are not included in this table (see SM2)."
    ), prop = note_prop)
  ))

out_path <- file.path(output_dir, "Scale_Reliability_Tables.docx")
print(doc, target = out_path)
message(sprintf("Saved: %s", basename(out_path)))
