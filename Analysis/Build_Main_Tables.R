# ============================================================
# Build_Main_Tables.R
# Generates APA7 Tables 1-3 for the BCAT manuscript main text.
#
# Input:  CSV files from the Results directory (Results.zip)
# Output: BCAT_Main_Tables.docx
#
# Table 1: Study Overview and Sample Characteristics
# Table 2: Arousal Gating by Awareness (H4B, Miss BF, H4C)
# Table 3: MAIA Validity (Confidence vs. Threshold)
# ============================================================

# Set Up ---------

# ── Config ────────────────────────────────────────────────────────────────────
FONT     <- "Times New Roman"
SZ       <- 12
SZ_NOTE  <- 10
RESULTS_DIR <- file.path(BASE_DIR, "Results")
OUTPUT_DIR  <- file.path(BASE_DIR, "Tables")

# ── Helpers ───────────────────────────────────────────────────────────────────

fmt_r <- function(x, d = 3) {
  # Format correlation: remove leading zero, preserve sign
  s <- formatC(abs(x), digits = d, format = "f")
  paste0(ifelse(x < 0, "-", ""), sub("^0", "", s))
}

fmt_p <- function(p) {
  ifelse(p < .001, "< .001", sub("^0", "", formatC(p, digits = 3, format = "f")))
}

fmt_bf <- function(bf) {
  ifelse(bf >= 10, formatC(bf, digits = 1, format = "f"),
         formatC(bf, digits = 2, format = "f"))
}

r_ci_fisher <- function(r, df) {
  # 95% CI via Fisher z; df = residual df (approx n_trials - 2)
  se <- 1 / sqrt(df - 1)
  z  <- atanh(r)
  c(lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se))
}

fmt_ci <- function(lo, hi, d = 3) {
  paste0("[", fmt_r(lo, d), ", ", fmt_r(hi, d), "]")
}

# APA7 three-rule flextable style
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

# Add a thin separator line before a specified row index
add_group_separator <- function(ft, i) {
  hline(ft, i = i - 1,
        border = officer::fp_border(width = 0.75, color = "black"),
        part = "body")
}

# Note paragraph appended after a table
table_note_par <- function(note_text, footnotes = NULL) {
  note_prop <- officer::fp_text(font.family = FONT, font.size = SZ_NOTE)
  bold_prop  <- officer::fp_text(font.family = FONT, font.size = SZ_NOTE, bold = TRUE)

  parts <- list(officer::ftext("Note. ", prop = bold_prop),
                officer::ftext(note_text, prop = note_prop))

  if (!is.null(footnotes)) {
    for (fn in footnotes) {
      parts <- c(parts,
                 list(officer::ftext(paste0(" ", fn$letter, " "), prop = bold_prop),
                      officer::ftext(fn$text, prop = note_prop)))
    }
  }
  do.call(officer::fpar, parts)
}

# ── Load data ─────────────────────────────────────────────────────────────────
aro   <- readr::read_csv(file.path(RESULTS_DIR, "table_arousal.csv"))
maia  <- readr::read_csv(file.path(RESULTS_DIR, "table_maia.csv"))
mbf   <- readr::read_csv(file.path(RESULTS_DIR, "miss_baseline_bf.csv"))
m4b_p <- readr::read_csv(file.path(RESULTS_DIR, "meta_h4b_pooled.csv"))
mmaia <- readr::read_csv(file.path(RESULTS_DIR, "meta_h3_maia_dissociation.csv"))

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 1: Study Overview and Sample Characteristics
# Hardcoded — values from hbd_summary.csv and Methods section.
# ══════════════════════════════════════════════════════════════════════════════

make_table_overview <- function() {
  yes <- "\u2713"
  no  <- "\u2014"

  d <- tibble::tibble(
    Feature = c(
      "N",
      "Age, M (SD)",
      "Women / Men / Other",
      "Setting",
      "Recruitment",
      "Staircase",
      "Arousal scale",
      "Salience manipulation",
      "Visual control (H4C)",
      "Physiological recording",
      "Preregistered"
    ),
    Study1A = c("181", "25.7 (8.4)",  "60 / 113 / 8",   "Online", "Prolific",   "Robbins\u2013Monro", "VAS (0\u201350)", no,  no,  no,  no),
    Study2  = c("166", "28.1 (10.3)", "51 / 114 / 1",   "Online", "Prolific",   "Robbins\u2013Monro", "VAS (0\u201350)", no,  no,  no,  no),
    Study3  = c("103", "19.5 (2.0)",  "76 / 27 / 0",    "Lab",    "UTM SONA",   "Fixed magnitudes",    "SAM (1\u20139)", yes, no,  no,  no),
    Study4  = c("131", "19.4 (3.4)",  "92 / 36 / 3",    "Online", "UTM SONA",   "Quest",               "Emoji (1\u20136)", yes, yes, no,  no),
    Study5  = c("206", "19.1 (2.1)",  "162 / 35 / 9",   "Lab",    "UTM SONA",   "Quest",               "Emoji (1\u20136)", yes, yes, yes, yes)
  )

  ft <- flextable::flextable(d) |>
    flextable::set_header_labels(
      Feature = "",
      Study1A = "Study 1A/1B\u1d43",
      Study2  = "Study 2",
      Study3  = "Study 3",
      Study4  = "Study 4",
      Study5  = "Study 5"
    ) |>
    flextable::bold(i = 1, part = "body") |>
    flextable::hline(
      i      = 3,
      border = officer::fp_border(width = 0.75, color = "black"),
      part   = "body"
    ) |>
    flextable::width(j = 1, width = 1.60) |>
    flextable::width(j = 2, width = 0.98) |>
    flextable::width(j = 3:6, width = 0.88) |>
    apa7_style()

  ft
}

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 2: Arousal Gating by Awareness
# ══════════════════════════════════════════════════════════════════════════════

build_table2 <- function() {
  study_labels <- c(Study1A = "1A", Study2 = "2", Study3 = "3",
                    Study4 = "4", Study5 = "5")
  bf_map <- setNames(mbf$BF01_multilevel, mbf$study)

  rows_study <- aro |>
    dplyr::mutate(
      study_label = study_labels[study],
      ci_h4b = purrr::map2(H4B_partial_r, H4B_df + 2, ~ r_ci_fisher(.x, .y)),
      ci_lo  = purrr::map_dbl(ci_h4b, "lo"),
      ci_hi  = purrr::map_dbl(ci_h4b, "hi"),
      ci_h4c = purrr::map2(
        ifelse(is.na(H4C_partial_r), 0, H4C_partial_r),
        ifelse(is.na(H4C_df), 100, H4C_df) + 2,
        ~ if (!is.na(.x) && .x != 0) r_ci_fisher(.x, .y) else list(lo = NA, hi = NA)
      ),
      ci_h4c_lo = purrr::map_dbl(ci_h4c, "lo"),
      ci_h4c_hi = purrr::map_dbl(ci_h4c, "hi"),
      bf_val     = bf_map[study]
    ) |>
    dplyr::transmute(
      Study    = paste("Study", study_label),
      n        = as.character(n),
      h4b_r    = fmt_r(H4B_partial_r),
      h4b_ci   = fmt_ci(ci_lo, ci_hi),
      bf       = dplyr::case_when(
        !is.na(bf_val) ~ fmt_bf(bf_val),
        TRUE ~ "\u2014"
      ),
      h4c_r    = dplyr::case_when(
        !is.na(H4C_partial_r) ~ fmt_r(H4C_partial_r),
        TRUE ~ "\u2014"
      ),
      h4c_ci   = dplyr::case_when(
        !is.na(ci_h4c_lo) ~ fmt_ci(ci_h4c_lo, ci_h4c_hi),
        TRUE ~ ""
      )
    )

  # Meta-analytic pooled row
  pm <- m4b_p[1, ]
  pooled_row <- tibble::tibble(
    Study  = "Pooled (meta-analytic)",
    n      = "k = 5",
    h4b_r  = fmt_r(pm$r_pooled),
    h4b_ci = fmt_ci(pm$r_lower, pm$r_upper),
    bf     = "\u2014",
    h4c_r  = "\u2014",
    h4c_ci = ""
  )

  list(
    data = dplyr::bind_rows(rows_study, pooled_row),
    n_study = nrow(rows_study)
  )
}

make_table2 <- function(obj) {
  df       <- obj$data
  n_study  <- obj$n_study

  flextable(df) |>
    set_header_labels(
      Study  = "Study", n = "n",
      h4b_r  = "partial r", h4b_ci = "95% CI",
      bf     = "BF\u2081\u2080 (null)\u1d43",
      h4c_r  = "partial r", h4c_ci = "95% CI"
    ) |>
    add_header_row(
      values     = c("", "", "Change \u00d7 Accuracy", "Miss vs.\nBaseline", "Breath > Visual"),
      colwidths  = c(1, 1, 2, 1, 2),
      top        = TRUE
    ) |>
    italic(i = 1, part = "header") |>
    width(j = "Study",  width = 1.8) |>
    width(j = "n",      width = 0.55) |>
    width(j = "h4b_r",  width = 0.85) |>
    width(j = "h4b_ci", width = 1.4) |>
    width(j = "bf",     width = 0.9) |>
    width(j = "h4c_r",  width = 0.85) |>
    width(j = "h4c_ci", width = 1.4) |>
    apa7_style() |>
    add_group_separator(n_study + 1) |>   # thin line before pooled row
    italic(i = nrow(df), part = "body")   # italicise pooled row label
}

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 2: MAIA Validity
# ══════════════════════════════════════════════════════════════════════════════

build_table3 <- function() {
  study_labels <- c(Study1A = "1A", Study2 = "2", Study3 = "3",
                    Study4 = "4", Study5 = "5")

  rows_study <- maia |>
    dplyr::mutate(study_label = study_labels[study]) |>
    dplyr::transmute(
      Study  = paste("Study", study_label),
      n      = as.character(n),
      h3a_r  = fmt_r(H3A_r),
      h3a_ci = fmt_ci(H3A_CI_lower, H3A_CI_upper),
      h3b_r  = fmt_r(H3B_r),
      h3b_ci = fmt_ci(H3B_CI_lower, H3B_CI_upper),
      bf     = dplyr::case_when(
        study == "Study3" ~ "\u2014\u1d47",   # footnote: no staircase threshold
        TRUE ~ formatC(H3B_BF01, digits = 2, format = "f")
      )
    )

  # Pooled rows from meta-analytic estimates
  conf_row <- mmaia |> dplyr::filter(grepl("Confidence", Contrast))
  thr_row  <- mmaia |> dplyr::filter(grepl("Threshold",  Contrast))

  pooled_row <- tibble::tibble(
    Study  = "Pooled (meta-analytic)",
    n      = paste0("k = ", c(conf_row$k_studies, thr_row$k_studies)[1]),
    h3a_r  = fmt_r(conf_row$r_pooled),
    h3a_ci = fmt_ci(conf_row$r_lower, conf_row$r_upper),
    h3b_r  = fmt_r(thr_row$r_pooled),
    h3b_ci = fmt_ci(thr_row$r_lower, thr_row$r_upper),
    bf     = ""
  )
  # Fix k for threshold (k=4)
  pooled_row$n <- paste0("k = ", conf_row$k_studies,
                         " / k = ", thr_row$k_studies)

  list(
    data    = dplyr::bind_rows(rows_study, pooled_row),
    n_study = nrow(rows_study)
  )
}

make_table3 <- function(obj) {
  df      <- obj$data
  n_study <- obj$n_study

  flextable(df) |>
    set_header_labels(
      Study  = "Study", n = "n",
      h3a_r  = "r", h3a_ci = "95% CI",
      h3b_r  = "r", h3b_ci = "95% CI",
      bf     = "BF\u2081\u2080 (null)"
    ) |>
    add_header_row(
      values    = c("", "", "MAIA \u2192 Confidence", "MAIA \u2192 Threshold", ""),
      colwidths = c(1, 1, 2, 2, 1),
      top       = TRUE
    ) |>
    italic(i = 1, part = "header") |>
    width(j = "Study",  width = 1.8) |>
    width(j = "n",      width = 0.65) |>
    width(j = "h3a_r",  width = 0.65) |>
    width(j = "h3a_ci", width = 1.55) |>
    width(j = "h3b_r",  width = 0.65) |>
    width(j = "h3b_ci", width = 1.55) |>
    width(j = "bf",     width = 0.95) |>
    apa7_style() |>
    add_group_separator(n_study + 1) |>
    italic(i = nrow(df), part = "body")
}

# ── Build Word document ───────────────────────────────────────────────────────
message("Building main tables Word document...")

ps_portrait <- officer::prop_section(
  page_size    = officer::page_size(width = 8.5, height = 11, orient = "portrait"),
  page_margins = officer::page_mar(top = 1, bottom = 1, left = 1, right = 1)
)

note_prop  <- officer::fp_text(font.family = FONT, font.size = SZ_NOTE)
bold_prop  <- officer::fp_text(font.family = FONT, font.size = SZ_NOTE, bold = TRUE)

add_table_block <- function(doc, table_num, title, ft, note_text, footnotes = NULL) {
  doc |>
    officer::body_add_par(sprintf("Table %d", table_num),
                          style = "heading 3") |>
    officer::body_add_par(title, style = "Normal") |>
    flextable::body_add_flextable(ft) |>
    officer::body_add_fpar(table_note_par(note_text, footnotes)) |>
    officer::body_add_par("", style = "Normal")
}

ft_overview <- make_table_overview()
obj2 <- build_table2()   # Arousal Gating -> Table 2
obj3 <- build_table3()   # MAIA Validity  -> Table 3

overview_note <- paste0(
  "UTM = University of Toronto Mississauga. ",
  "VAS = visual analogue scale (0\u201350). ",
  "SAM = Self-Assessment Manikin (1\u20139). ",
  "Quest = Watson\u2013Pelli adaptive Bayesian staircase (75% correct target). ",
  "Robbins\u2013Monro = stochastic approximation staircase. ",
  "Salience manipulation = high (abrupt) vs. low (gradual) onset of change. ",
  "Visual control = between-groups visual pacing condition (interoceptive specificity, H4C). ",
  "\u1d43Study 1B was conducted in the same session as Study 1A using an ascending-limits ",
  "procedure (N = 181, same sample). Study 1B data are used for procedure comparison only ",
  "(Supplementary S1.1) and excluded from all arousal and MAIA analyses."
)

doc <- officer::read_docx() |>
  officer::body_set_default_section(ps_portrait) |>

  # ── Table 1 (Study Overview) ──────────────────────────────────────────────
  officer::body_add_par("Table 1", style = "heading 3") |>
  officer::body_add_par(
    "Study Overview and Sample Characteristics",
    style = "Normal") |>
  flextable::body_add_flextable(ft_overview) |>
  officer::body_add_fpar(table_note_par(overview_note)) |>
  officer::body_add_par("", style = "Normal") |>

  # ── Table 2 (Arousal Gating) ──────────────────────────────────────────────
  add_table_block(
    2,
    "Arousal Gating by Awareness: Change \u00d7 Accuracy Interaction Across Studies",
    make_table2(obj2),
    paste0(
      "Change \u00d7 Accuracy partial r reflects the moderation of the breathing rate ",
      "change\u2013arousal relationship by trial-level detection accuracy (Hit = detected, ",
      "Miss = not detected), from linear mixed-effects models with random intercepts per participant. ",
      "Negative values indicate stronger change\u2013arousal coupling on detected than missed trials. ",
      "Likelihood ratio tests confirmed the quadratic Change\u00b2 term significantly predicted ",
      "subjective arousal in all five samples (all p < .001). ",
      "Breath > Visual partial r indexes the arousal advantage for the Breath over Visual condition; ",
      "available in Studies 4 and 5 only. ",
      "Pooled row reports a random-effects meta-analytic estimate (k = 5; I\u00b2 = 93.8\u0025; ",
      "heterogeneity driven by arousal scale differences across studies). ",
      "Dashes indicate the measure was not available for that study."
    ),
    footnotes = list(
      list(
        letter = "a",
        text = paste0(
          "BF\u2081\u2080 (null) = Bayes factor favouring the null that missed-change trials produce no ",
          "arousal above the no-change baseline, from Bayesian multilevel models with default JZS prior ",
          "(r = .707). Study 3 had no no-change baseline condition and is excluded."
        )
      )
    )
  ) |>

  # ── Table 3 (MAIA Validity) ───────────────────────────────────────────────
  add_table_block(
    3,
    "MAIA Validity: Dissociation of Confidence Prediction from Detection Threshold Prediction",
    make_table3(obj3),
    paste0(
      "MAIA \u2192 Confidence: Pearson r between MAIA total score and mean trial-level detection ",
      "confidence. MAIA \u2192 Threshold: Pearson r between MAIA total score and individually ",
      "estimated detection threshold (staircase-derived). ",
      "BF\u2081\u2080 (null) = Bayes factor for the null MAIA\u2013threshold association, ",
      "default JZS prior (r = .707). ",
      "Pooled rows report random-effects meta-analytic estimates ",
      "(k = 5 for MAIA\u2013confidence; k = 4 for MAIA\u2013threshold). ",
      "95% CIs based on Fisher r-to-z transformation."
    ),
    footnotes = list(
      list(
        letter = "a",
        text = paste0(
          "Study 3 used fixed change magnitudes; mean detection accuracy served as the ",
          "interoceptive sensitivity proxy rather than a staircase-derived threshold. ",
          "BF\u2081\u2080 = 0.85 (inconclusive). Excluded from the pooled threshold estimate."
        )
      )
    )
  )

out_path <- file.path(OUTPUT_DIR, "BCAT_Main_Tables.docx")
print(doc, target = out_path)
message(sprintf("Saved: %s", basename(out_path)))
