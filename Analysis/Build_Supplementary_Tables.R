# ============================================================
# Build_Supplementary_Tables.R
# Generates supplementary results tables for the BCAT manuscript.
#
# Input:  CSV files from the Results directory
# Output: BCAT_Supplementary_Tables.docx
#
# ST1: Detection Sensitivity -- Full Logistic Model Results
# ST3: Threshold Validation (H1: Direction, H2: Salience, H5: Predictive Validity)
# ST4: Staircase Convergence
# ST5: Threshold Descriptives by Salience x Direction
# ST6: Test Block d' (Studies 4 and 5)
# ST7: Test-Retest Reliability (Study 5 ICC)
# ST8: Study 5 Belt Compliance and Physiological Arousal
# ST8: MAIA-Confidence Partial Correlations Controlling for Self-Esteem (Study 5)
# ============================================================

# Set Up ---------

# ── Config ────────────────────────────────────────────────────────────────────
FONT      <- "Times New Roman"
SZ        <- 11
SZ_NOTE   <- 9
RESULTS_DIR <- file.path(BASE_DIR, "Results")
OUTPUT_DIR  <- file.path(BASE_DIR, "Tables")

# ── Helpers ───────────────────────────────────────────────────────────────────
fmt_r3 <- function(x, d = 3) {
  s <- formatC(abs(x), digits = d, format = "f")
  paste0(ifelse(x < 0, "-", ""), sub("^0", "", s))
}

fmt_r2 <- function(x) fmt_r3(x, d = 2)

fmt_p <- function(p) {
  ifelse(p < .001, "< .001", sub("^0", "", formatC(p, digits = 3, format = "f")))
}

fmt_ci <- function(lo, hi, d = 3) {
  paste0("[", fmt_r3(lo, d), ", ", fmt_r3(hi, d), "]")
}

r_ci <- function(r, df) {
  se <- 1 / sqrt(df - 1)
  z  <- atanh(r)
  c(lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se))
}

# APA7 three-rule style
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

thin_rule <- function(ft, i) {
  hline(ft, i = i, border = officer::fp_border(width = 0.75, color = "black"), part = "body")
}

note_par <- function(note_text, footnotes = NULL) {
  np  <- officer::fp_text(font.family = FONT, font.size = SZ_NOTE)
  bld <- officer::fp_text(font.family = FONT, font.size = SZ_NOTE, bold = TRUE)
  parts <- list(officer::ftext("Note. ", prop = bld),
                officer::ftext(note_text, prop = np))
  if (!is.null(footnotes)) {
    for (fn in footnotes) {
      parts <- c(parts,
                 list(officer::ftext(paste0(" ", fn$l, " "), prop = bld),
                      officer::ftext(fn$t, prop = np)))
    }
  }
  do.call(officer::fpar, parts)
}

add_table_block <- function(doc, label, title, ft, note, footnotes = NULL) {
  doc |>
    officer::body_add_par(label, style = "heading 3") |>
    officer::body_add_par(title, style = "Normal") |>
    flextable::body_add_flextable(ft) |>
    officer::body_add_fpar(note_par(note, footnotes)) |>
    officer::body_add_par("", style = "Normal")
}

# fmt_b_se: "b (SE)" formatted string; "--" for NA
fmt_b_se <- function(b, se, d = 2) {
  ifelse(is.na(b) | is.na(se), "\u2014",
         sprintf("%s (%s)",
                 formatC(round(b, d), digits = d, format = "f"),
                 formatC(round(se, d), digits = d, format = "f")))
}

# ── Load data ─────────────────────────────────────────────────────────────────
det   <- readr::read_csv(file.path(RESULTS_DIR, "table_detection_change2.csv"))
aro   <- readr::read_csv(file.path(RESULTS_DIR, "table_arousal.csv"))
conv  <- readr::read_csv(file.path(RESULTS_DIR, "table_staircase_convergence.csv"))
thr   <- readr::read_csv(file.path(RESULTS_DIR, "table_threshold_descriptives.csv"))
dprime    <- readr::read_csv(file.path(RESULTS_DIR, "test_block_dprime_by_group.csv"))
dprime_pc <- readr::read_csv(file.path(RESULTS_DIR, "test_block_accuracy_descriptives.csv"))
rel   <- readr::read_csv(file.path(RESULTS_DIR, "table_reliability.csv"))
belt  <- readr::read_csv(file.path(RESULTS_DIR, "table_belt_physio_arousal.csv"))
pcor  <- readr::read_csv(file.path(RESULTS_DIR, "s7_partial_correlations.csv"))
mlm   <- readr::read_csv(file.path(RESULTS_DIR, "s7_multilevel_confidence_results.csv"))
id4   <- readr::read_csv(file.path(RESULTS_DIR, "id_correlations_study4.csv"))
id5   <- readr::read_csv(file.path(RESULTS_DIR, "id_correlations_study5.csv"))
val   <- readr::read_csv(file.path(RESULTS_DIR, "table_validation.csv"))

# Residual df proxy for partial r in logistic detection models
df_proxy <- setNames(aro$H4B_df, aro$study)

# ══════════════════════════════════════════════════════════════════════════════
# ST1: Detection Sensitivity -- Full Logistic Model Results
# ══════════════════════════════════════════════════════════════════════════════

build_st1 <- function() {
  study_map <- c(Study1_TaskA = "Study1A", Study2 = "Study2",
                 Study3 = "Study3", Study4 = "Study4", Study5 = "Study5")
  label_map <- c(Study1_TaskA = "1A (staircase)", Study2 = "2",
                 Study3 = "3 (fixed magnitudes)", Study4 = "4", Study5 = "5")
  
  det |>
    dplyr::filter(study != "Study1_TaskB") |>
    dplyr::mutate(
      df_resid  = df_proxy[study_map[study]],
      partial_r = z_Change2 / sqrt(z_Change2^2 + df_resid),
      ci_lo     = purrr::map2_dbl(partial_r, df_resid + 2, ~ r_ci(.x, .y)["lo"]),
      ci_hi     = purrr::map2_dbl(partial_r, df_resid + 2, ~ r_ci(.x, .y)["hi"])
    ) |>
    dplyr::filter(!is.na(df_resid)) |>
    dplyr::transmute(
      Study     = paste("Study", label_map[study]),
      n         = as.character(n),
      OR_Change  = formatC(OR_Change,  digits = 2, format = "f"),
      OR_Change2 = formatC(OR_Change2, digits = 1, format = "f"),
      partial_r  = fmt_r3(partial_r),
      ci         = fmt_ci(ci_lo, ci_hi),
      lrt_p      = fmt_p(lrt_p_quad)
    )
}

make_st1 <- function(df) {
  flextable(df) |>
    set_header_labels(
      Study     = "Study",
      n         = "n",
      OR_Change  = "OR (Change)",
      OR_Change2 = "OR (Change\u00b2)",
      partial_r  = "partial r",
      ci         = "95% CI",
      lrt_p      = "p (LRT)"
    ) |>
    add_header_row(
      values    = c("", "", "Odds Ratios", "Quadratic Term (Change\u00b2)", ""),
      colwidths = c(1, 1, 2, 2, 1),
      top       = TRUE
    ) |>
    italic(i = 1, part = "header") |>
    width(j = "Study",     width = 1.8) |>
    width(j = "n",         width = 0.5) |>
    width(j = "OR_Change",  width = 1.1) |>
    width(j = "OR_Change2", width = 1.1) |>
    width(j = "partial_r", width = 0.85) |>
    width(j = "ci",        width = 1.5) |>
    width(j = "lrt_p",     width = 0.85) |>
    apa7_style()
}

# ══════════════════════════════════════════════════════════════════════════════
# ST2: Threshold Validation (H1: Direction, H2: Salience, H5: Predictive Validity)
# Source: table_validation.csv
# ══════════════════════════════════════════════════════════════════════════════

build_st2_val <- function() {
  study_labels <- c(Study1A = "1A", Study1B = "1B\u1d43",
                    Study4 = "4", Study5 = "5")
  val |>
    dplyr::filter(study %in% names(study_labels)) |>
    dplyr::mutate(
      Study  = paste("Study", study_labels[study]),
      h1_bse = fmt_b_se(H1_b, H1_SE),
      h1_r   = dplyr::if_else(is.na(H1_partial_r), "\u2014", fmt_r2(H1_partial_r)),
      h1_p   = dplyr::if_else(is.na(H1_p),         "\u2014", fmt_p(H1_p)),
      h2_bse = dplyr::if_else(is.na(H2_b),         "\u2014", fmt_b_se(H2_b, H2_SE)),
      h2_r   = dplyr::if_else(is.na(H2_partial_r), "\u2014", fmt_r2(H2_partial_r)),
      h2_p   = dplyr::if_else(is.na(H2_p),         "\u2014", fmt_p(H2_p)),
      h5_bse = dplyr::if_else(is.na(H5_b),         "\u2014", fmt_b_se(H5_b, H5_SE)),
      h5_p   = dplyr::if_else(is.na(H5_p),         "\u2014", fmt_p(H5_p))
    ) |>
    dplyr::transmute(Study, h1_bse, h1_r, h1_p, h2_bse, h2_r, h2_p, h5_bse, h5_p)
}

make_st2_val <- function(df) {
  flextable(df) |>
    set_header_labels(
      Study  = "Study",
      h1_bse = "b (SE)", h1_r = "partial r", h1_p = "p",
      h2_bse = "b (SE)", h2_r = "partial r", h2_p = "p",
      h5_bse = "b (SE)", h5_p = "p"
    ) |>
    add_header_row(
      values    = c("", "H1: Direction", "H2: Salience", "H5: Test validity"),
      colwidths = c(1, 3, 3, 2),
      top       = TRUE
    ) |>
    italic(i = 1, part = "header") |>
    width(j = "Study",                          width = 0.80) |>
    width(j = c("h1_bse", "h2_bse", "h5_bse"), width = 1.00) |>
    width(j = c("h1_r", "h2_r"),               width = 0.70) |>
    width(j = c("h1_p", "h2_p", "h5_p"),       width = 0.55) |>
    apa7_style()
}

# ══════════════════════════════════════════════════════════════════════════════
# ST3: Staircase Convergence
# ══════════════════════════════════════════════════════════════════════════════

build_st3 <- function() {
  conv |>
    dplyr::mutate(
      pct_reduction = formatC(pct_reduction, digits = 1, format = "f"),
      cohens_dz     = fmt_r3(cohens_dz),
      p_str         = fmt_p(p),
      cell_label    = dplyr::case_when(
        cell == "Overall"    ~ "All conditions",
        TRUE ~ gsub("x", " \u00d7 ", cell)
      )
    ) |>
    dplyr::transmute(
      Study    = study,
      Cell     = cell_label,
      n        = as.character(n),
      M_trials = formatC(n_trials_mean, digits = 1, format = "f"),
      SD_first6 = formatC(SD_first6, digits = 3, format = "f"),
      SD_last6  = formatC(SD_last6,  digits = 3, format = "f"),
      pct_red  = paste0(pct_reduction, "%"),
      d_z      = cohens_dz,
      p        = p_str
    )
}

make_st3 <- function(df) {
  # Find rows where study changes for thin separators
  study_breaks <- which(df$Study != dplyr::lag(df$Study, default = ""))[-1] - 1
  
  ft <- flextable(df) |>
    set_header_labels(
      Study    = "Study",
      Cell     = "Condition",
      n        = "n",
      M_trials = "Mean trials",
      SD_first6 = "SD (first 6)",
      SD_last6  = "SD (last 6)",
      pct_red  = "Reduction",
      d_z      = "d\u209a",
      p        = "p"
    ) |>
    merge_v(j = "Study") |>
    valign(j = "Study", valign = "top", part = "body") |>
    width(j = "Study",    width = 0.7) |>
    width(j = "Cell",     width = 1.6) |>
    width(j = "n",        width = 0.45) |>
    width(j = "M_trials", width = 0.85) |>
    width(j = c("SD_first6","SD_last6"), width = 0.95) |>
    width(j = "pct_red",  width = 0.85) |>
    width(j = "d_z",      width = 0.7) |>
    width(j = "p",        width = 0.7) |>
    apa7_style()
  
  for (i in study_breaks) ft <- thin_rule(ft, i)
  ft
}

# ══════════════════════════════════════════════════════════════════════════════
# ST4: Threshold Descriptives by Salience x Direction
# ══════════════════════════════════════════════════════════════════════════════

build_st4 <- function() {
  thr |>
    dplyr::mutate(
      ses_label = dplyr::case_when(
        is.na(ses) | ses == "ses1" ~ "Session 1",
        ses == "ses2"              ~ "Session 2",
        TRUE                       ~ ses
      )
    ) |>
    dplyr::transmute(
      Study     = paste(study, ses_label),
      Salience  = Salience,
      Direction = Direction,
      n         = as.character(n_ppts),
      M         = formatC(M_threshold * 100, digits = 1, format = "f"),
      SD        = formatC(SD_threshold * 100, digits = 1, format = "f")
    )
}

make_st4 <- function(df) {
  breaks <- which(df$Study != dplyr::lag(df$Study, default = ""))[-1] - 1
  
  ft <- flextable(df) |>
    set_header_labels(
      Study = "Study / Session", Salience = "Salience",
      Direction = "Direction", n = "n",
      M = "M (%)", SD = "SD (%)"
    ) |>
    merge_v(j = c("Study","Salience")) |>
    valign(j = c("Study","Salience"), valign = "top", part = "body") |>
    width(j = "Study",     width = 1.6) |>
    width(j = "Salience",  width = 0.8) |>
    width(j = "Direction", width = 0.85) |>
    width(j = "n",         width = 0.5) |>
    width(j = "M",         width = 0.7) |>
    width(j = "SD",        width = 0.7) |>
    apa7_style()
  
  for (i in breaks) ft <- thin_rule(ft, i)
  ft
}

# ══════════════════════════════════════════════════════════════════════════════
# ST5: Test Block d' (Studies 4 and 5) — direction-stratified
# ══════════════════════════════════════════════════════════════════════════════

build_st5 <- function() {
  cond_labels <- c(
    S4_Breath       = "Study 4 / Breath",
    S4_Visual       = "Study 4 / Visual",
    S5_ses1_Breath  = "Study 5 / ses1 Breath",
    S5_ses1_Visual  = "Study 5 / ses1 Visual",
    S5_ses2_Breath  = "Study 5 / ses2 Breath"
  )
  
  # Coerce join keys to character in both tables to avoid factor mismatch.
  # Drop n_trials/n_pp/Pc from dprime first so the join works whether or not
  # test_block_accuracy.R has already embedded them in the CSV.
  dprime_clean <- dprime |>
    dplyr::mutate(dplyr::across(c(ses_cond_f, salience_f, direction_f), as.character)) |>
    dplyr::select(-dplyr::any_of(c("n_trials", "n_pp", "Pc")))
  
  pc_clean <- dprime_pc |>
    dplyr::mutate(dplyr::across(c(ses_cond_f, salience_f, direction_f), as.character)) |>
    dplyr::select(ses_cond_f, salience_f, direction_f, n_trials, n_pp, Pc)
  
  dprime_clean |>
    dplyr::left_join(pc_clean, by = c("ses_cond_f", "salience_f", "direction_f")) |>
    dplyr::mutate(
      Condition = dplyr::recode(ses_cond_f, !!!cond_labels),
      Condition = factor(Condition, levels = unname(cond_labels))
    ) |>
    dplyr::arrange(Condition, salience_f, direction_f) |>
    dplyr::transmute(
      Condition = as.character(Condition),
      Salience  = salience_f,
      Direction = direction_f,
      N_trials  = as.character(n_trials),
      N_ppts    = as.character(n),
      Pc        = formatC(Pc,          digits = 3, format = "f"),
      dprime    = formatC(mean_dprime, digits = 2, format = "f"),
      dprime_sd = formatC(sd_dprime,   digits = 2, format = "f")
    )
}

make_st5 <- function(df) {
  # Thin rule between conditions
  cond_breaks <- which(df$Condition != dplyr::lag(df$Condition,
                                                  default = ""))[-1] - 1
  
  ft <- flextable(df) |>
    set_header_labels(
      Condition = "Condition",
      Salience  = "Salience",
      Direction = "Direction",
      N_trials  = "N trials",
      N_ppts    = "N participants",
      Pc        = "P(correct)",
      dprime    = "d\u2019",
      dprime_sd = "SD(d\u2019)"
    ) |>
    merge_v(j = c("Condition", "Salience")) |>
    valign(j = c("Condition", "Salience"), valign = "top", part = "body") |>
    width(j = "Condition",  width = 1.35) |>
    width(j = "Salience",   width = 0.70) |>
    width(j = "Direction",  width = 0.75) |>
    width(j = "N_trials",   width = 0.65) |>
    width(j = "N_ppts",     width = 1.00) |>
    width(j = "Pc",         width = 0.80) |>
    width(j = "dprime",     width = 0.55) |>
    width(j = "dprime_sd",  width = 0.65) |>
    apa7_style()
  
  for (i in cond_breaks) ft <- thin_rule(ft, i)
  ft
}

# ══════════════════════════════════════════════════════════════════════════════
# ST6: Test-Retest Reliability (Study 5 ICC)
# ══════════════════════════════════════════════════════════════════════════════

build_st6 <- function() {
  rel |>
    dplyr::mutate(
      cond_label  = gsub("_", " / ", condition),
      group_short = dplyr::case_when(
        grepl("Breath-first", group) ~ "Breath\u2013Breath",
        grepl("Visual-first", group) ~ "Visual\u2013Breath",
        TRUE ~ group
      )
    ) |>
    dplyr::transmute(
      Group     = group_short,
      Condition = cond_label,
      n         = as.character(n),
      ICC       = fmt_r3(icc),
      CI        = fmt_ci(icc_lower, icc_upper)
    )
}

make_st6 <- function(df) {
  breaks <- which(df$Group != dplyr::lag(df$Group, default = ""))[-1] - 1
  
  ft <- flextable(df) |>
    set_header_labels(
      Group = "Group", Condition = "Condition",
      n = "n", ICC = "ICC(2,1)", CI = "95% CI"
    ) |>
    merge_v(j = "Group") |>
    valign(j = "Group", valign = "top", part = "body") |>
    width(j = "Group",     width = 1.4) |>
    width(j = "Condition", width = 1.5) |>
    width(j = "n",         width = 0.45) |>
    width(j = "ICC",       width = 0.7) |>
    width(j = "CI",        width = 1.5) |>
    apa7_style()
  
  for (i in breaks) ft <- thin_rule(ft, i)
  ft
}

# ══════════════════════════════════════════════════════════════════════════════
# ST7: Belt Compliance and Physiological Arousal (Study 5)
# ══════════════════════════════════════════════════════════════════════════════

build_st7 <- function() {
  belt |>
    dplyr::transmute(
      Condition = condition,
      b         = formatC(estimate,   digits = 4, format = "f"),
      SE        = formatC(std.error,  digits = 4, format = "f"),
      t_or_chi  = dplyr::case_when(
        !is.na(statistic) ~ formatC(statistic, digits = 2, format = "f"),
        TRUE              ~ formatC(lrt_chi2,  digits = 2, format = "f")
      ),
      p         = fmt_p(p.value),
      CI_95     = dplyr::case_when(
        !is.na(conf.low) ~ fmt_ci(conf.low, conf.high, d = 4),
        TRUE             ~ ""
      ),
      note_col  = note
    ) |>
    dplyr::select(-note_col)
}

make_st7 <- function(df) {
  flextable(df) |>
    set_header_labels(
      Condition = "Model / Condition",
      b = "b", SE = "SE",
      t_or_chi = "t / \u03c7\u00b2",
      p = "p",
      CI_95 = "95% CI"
    ) |>
    width(j = "Condition", width = 2.8) |>
    width(j = "b",         width = 0.8) |>
    width(j = "SE",        width = 0.7) |>
    width(j = "t_or_chi",  width = 0.7) |>
    width(j = "p",         width = 0.7) |>
    width(j = "CI_95",     width = 1.6) |>
    apa7_style()
}

# ══════════════════════════════════════════════════════════════════════════════
# ST8: MAIA-Confidence Partial Correlations (Study 5, S8)
# ══════════════════════════════════════════════════════════════════════════════

build_st8 <- function() {
  # CSV first column is 'label' (values: "ses1", "ses2"); rename for clarity.
  # Each non-null control appears twice: once for single-control model and
  # once for the both-controls model (rows 3-4 share the same r value).
  # Use row_number() within (ses, control) to distinguish them.
  pcor |>
    dplyr::rename(ses = label) |>
    dplyr::filter(predictor == "MAIA_total_z") |>
    dplyr::group_by(ses, control) |>
    dplyr::mutate(occurrence = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      Session = dplyr::case_when(ses == "ses1" ~ "1", ses == "ses2" ~ "2"),
      Control = dplyr::case_when(
        control == "none"             ~ "None",
        control == "MSES_total_z"     & occurrence == 1L ~ "MSES total",
        control == "MSES_selfdoubt_z" & occurrence == 1L ~ "MSES self-doubt",
        occurrence == 2L              ~ "Both",
        TRUE                          ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(Control)) |>
    # occurrence == 2 produces two "Both" rows (same r); keep first
    dplyr::distinct(Session, Control, .keep_all = TRUE) |>
    dplyr::arrange(Session,
                   match(Control, c("None", "MSES total", "MSES self-doubt", "Both"))) |>
    dplyr::transmute(
      Session   = Session,
      Control   = Control,
      n         = as.character(n),
      partial_r = fmt_r3(r),
      p         = fmt_p(p)
    )
}

make_st8 <- function(df) {
  breaks <- which(df$Session != dplyr::lag(df$Session, default = ""))[-1] - 1
  
  ft <- flextable(df) |>
    set_header_labels(
      Session   = "Session",
      Control   = "Covariate controlled",
      n         = "n",
      partial_r = "partial r",
      p         = "p"
    ) |>
    merge_v(j = "Session") |>
    valign(j = "Session", valign = "top", part = "body") |>
    width(j = "Session",   width = 0.7) |>
    width(j = "Control",   width = 1.7) |>
    width(j = "n",         width = 0.5) |>
    width(j = "partial_r", width = 0.9) |>
    width(j = "p",         width = 0.7) |>
    apa7_style()
  
  for (i in breaks) ft <- thin_rule(ft, i)
  ft
}

# ── Build Word document ───────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
# ST9: Individual Differences -- Directionally Replicating Findings (Studies 4-5)
# Source: id_correlations_study4.csv, id_correlations_study5.csv
# Shows only the five pairs that replicated directionally across both studies.
# All correlations uncorrected for multiple comparisons; treat as exploratory.
# ══════════════════════════════════════════════════════════════════════════════

build_st9_id <- function() {
  # The five directionally-replicating task x questionnaire pairs
  target_pairs <- tibble::tibble(
    task_feature  = c("mean_Confidence", "mean_Confidence",
                      "slope_hit_overall", "conf_acc_r", "mean_arousal_overall"),
    questionnaire = c("MAIA_total", "MAIA_EmoAware",
                      "MAIA_EmoAware", "MAIA_Trusting", "MAIA_total")
  )
  
  feature_labels <- c(
    mean_Confidence    = "Mean detection confidence",
    slope_hit_overall  = "Arousal slope: hit trials",
    conf_acc_r         = "Confidence\u2013accuracy correlation",
    mean_arousal_overall = "Mean arousal"
  )
  quest_labels <- c(
    MAIA_total    = "MAIA Total",
    MAIA_EmoAware = "MAIA Emotional Awareness",
    MAIA_Trusting = "MAIA Trusting"
  )
  
  s4 <- id4 |>
    dplyr::semi_join(target_pairs, by = c("task_feature", "questionnaire")) |>
    dplyr::transmute(
      task_feature, questionnaire,
      n4  = as.character(n),
      r4  = fmt_r2(r),
      ci4 = fmt_ci(ci_lower, ci_upper, d = 2),
      p4  = fmt_p(p)
    )
  
  s5 <- id5 |>
    dplyr::semi_join(target_pairs, by = c("task_feature", "questionnaire")) |>
    dplyr::transmute(
      task_feature, questionnaire,
      n5  = as.character(n),
      r5  = fmt_r2(r),
      ci5 = fmt_ci(ci_lower, ci_upper, d = 2),
      p5  = fmt_p(p)
    )
  
  target_pairs |>
    dplyr::left_join(s4, by = c("task_feature", "questionnaire")) |>
    dplyr::left_join(s5, by = c("task_feature", "questionnaire")) |>
    dplyr::mutate(
      Task        = feature_labels[task_feature],
      Questionnaire = quest_labels[questionnaire]
    ) |>
    dplyr::select(Task, Questionnaire, n4, r4, ci4, p4, n5, r5, ci5, p5)
}

make_st9_id <- function(df) {
  flextable(df) |>
    set_header_labels(
      Task          = "Task feature",
      Questionnaire = "Questionnaire",
      n4 = "n", r4 = "r", ci4 = "95% CI", p4 = "p",
      n5 = "n", r5 = "r", ci5 = "95% CI", p5 = "p"
    ) |>
    add_header_row(
      values    = c("", "", "Study 4", "Study 5"),
      colwidths = c(1, 1, 4, 4),
      top       = TRUE
    ) |>
    italic(i = 1, part = "header") |>
    merge_v(j = "Task") |>
    valign(j = "Task", valign = "top", part = "body") |>
    width(j = "Task",          width = 1.90) |>
    width(j = "Questionnaire", width = 1.55) |>
    width(j = c("n4","n5"),    width = 0.35) |>
    width(j = c("r4","r5"),    width = 0.40) |>
    width(j = c("ci4","ci5"),  width = 1.10) |>
    width(j = c("p4","p5"),    width = 0.50) |>
    apa7_style()
}

message("Building supplementary tables Word document...")

ps <- officer::prop_section(
  page_size    = officer::page_size(width = 8.5, height = 11, orient = "portrait"),
  page_margins = officer::page_mar(top = 1, bottom = 1, left = 1, right = 1)
)

st1_df <- build_st1()
st2_val_df <- build_st2_val()
st3_df <- build_st3()
st4_df <- build_st4()
st5_df <- build_st5()
st6_df <- build_st6()
st7_df <- build_st7()
st8_df <- build_st8()
st9_id_df <- build_st9_id()

doc <- officer::read_docx() |>
  officer::body_set_default_section(ps) |>
  
  add_table_block(
    "Table ST1",
    "Detection Sensitivity: Logistic Mixed-Effects Model Results Across Studies",
    make_st1(st1_df),
    paste0(
      "Multilevel logistic regression predicting 3AFC detection accuracy from ",
      "linear (Change) and quadratic (Change\u00b2) breathing rate change terms, ",
      "with random intercepts per participant. ",
      "OR = odds ratio for a one-unit increase in Change or Change\u00b2. ",
      "partial r for Change\u00b2 is approximated from the z-statistic using the arousal model ",
      "residual df as a proxy for trial count. ",
      "p (LRT) = likelihood ratio test p-value for the quadratic term improvement. ",
      "Study 1B (ascending-limits procedure) is excluded as structurally incompatible."
    )
  ) |>
  
  add_table_block(
    "Table ST2",
    "Threshold Validation: Direction, Salience, and Predictive Validity",
    make_st2_val(st2_val_df),
    paste0(
      "H1: effect of breathing direction (Slower vs. Faster) on individually estimated staircase threshold, ",
      "from mixed-effects regression with random intercepts per participant. ",
      "H2: effect of salience manipulation (High vs. Low abruptness of onset) on threshold; Studies 4 and 5 only. ",
      "H5: logistic regression of staircase-derived threshold on test-block accuracy (predictive validity); ",
      "Studies 4 and 5 only. ",
      "b and SE are unstandardised coefficients. ",
      "partial r computed from the t-statistic (r = t / \u221a(t\u00b2 + df)). ",
      "Dashes indicate the hypothesis was not tested in that study."
    ),
    footnotes = list(list(
      l = "a",
      t = paste0("Study 1B used an ascending-limits procedure rather than a staircase. ",
                 "H1 is reported for comparison purposes; ascending-limits thresholds are not comparable ",
                 "with staircase thresholds from other studies.")
    ))
  ) |>
  
  add_table_block(
    "Table ST3",
    "Staircase Convergence: Within-Person Variability Reduction Across Studies",
    make_st3(st3_df),
    paste0(
      "Within-person SD of the staircase change magnitude was computed separately for ",
      "the first and last 6 effective change trials per staircase condition. ",
      "Reduction = percentage decrease from SD (first 6) to SD (last 6). ",
      "d\u209a = Cohen\u2019s d for paired differences (within-person). ",
      "All reductions significant at p < .001. ",
      "Study 2 used a single undifferentiated staircase (all conditions pooled); ",
      "Studies 4 and 5 report separate estimates per Salience \u00d7 Direction cell."
    )
  ) |>
  
  add_table_block(
    "Table ST4",
    "Threshold Descriptives by Salience and Direction (Studies 4 and 5)",
    make_st4(st4_df),
    paste0(
      "Detection thresholds expressed as percent change in breathing rate from baseline. ",
      "Thresholds estimated as the mean of the final two Quest staircase trials per condition. ",
      "Higher threshold values indicate lower detection sensitivity (larger change required for ",
      "reliable detection). Study 5 Session 2 = Block 2 (all participants completed Breath condition)."
    )
  ) |>
  
  add_table_block(
    "Table ST5",
    "Test Block Detection Sensitivity (d\u2019) by Condition, Salience, and Direction",
    make_st5(st5_df),
    paste0(
      "d\u2019 computed via numerical inversion of the unbiased 3AFC operating characteristic; ",
      "SD(d\u2019) = within-condition standard deviation across participants. ",
      "Test block thresholds were set separately per direction as the mean of the ",
      "high- and low-salience staircase estimates for that direction. ",
      "High-salience d\u2019 approximates the Quest target (\u22481.50 / 75% correct) in most conditions. ",
      "ses1 = Session 1; ses2 = Session 2 (Study 5 only; all participants completed Breath condition). ",
      "Study 4 Visual group was absent from prior versions of this table and has been restored."
    )
  ) |>
  
  add_table_block(
    "Table ST6",
    "Test-Retest Reliability of Staircase Threshold Estimates (Study 5)",
    make_st6(st6_df),
    paste0(
      "ICC(2,1) = two-way random-effects intraclass correlation, single rater, absolute agreement. ",
      "Breath\u2013Breath group: participants who completed Breath condition in both blocks (n \u2248 65). ",
      "Visual\u2013Breath group: participants who completed Visual condition in Block 1 and Breath in Block 2; ",
      "ICC reflects cross-condition (not test-retest) reliability. ",
      "Pre-registered criterion ICC \u2265 .70 was not met at the condition level. ",
      "Aggregate r across conditions: Breath\u2013Breath = .673, Visual\u2013Breath = .318 (see main text)."
    )
  ) |>
  
  add_table_block(
    "Table ST7",
    "Belt Compliance and Physiological Arousal Effects (Study 5)",
    make_st7(st7_df),
    paste0(
      "Rows 1\u20132: multilevel linear regression of within-person respiratory duration change ",
      "(dur_b1_vs_b4) on subjective arousal, estimated separately for detected (Hits) and ",
      "missed (Misses) trials. Negative b = larger breath duration change predicts higher arousal. ",
      "Row 3: likelihood ratio test (LRT) for the duration \u00d7 Accuracy interaction, testing whether ",
      "detection moderates the physiology\u2013arousal relationship. ",
      "Analysis restricted to participants with non-unusable belt quality and valid session onset ",
      "(N = 171; 9,270 paced trials)."
    )
  ) |>
  
  add_table_block(
    "Table ST8",
    "MAIA-Confidence Association Controlling for Self-Esteem (Study 5)",
    make_st8(st8_df),
    paste0(
      "Partial correlations of MAIA total score with session-level mean detection confidence, ",
      "after controlling for self-esteem measures (Study 5). ",
      "Session 1 = Block 1 (Breath or Visual condition, counterbalanced; N = 201). ",
      "Session 2 = Block 2 (all participants: Breath condition; N = 155). ",
      "MSES = Multidimensional Self-Esteem Scale\u201312 (Rentzsch et al., 2022). ",
      "The MAIA\u2013confidence association remained significant and largely unchanged in magnitude ",
      "after controlling for both MSES total and self-doubt subscale."
    )
  ) |>
  
  add_table_block(
    "Table ST9",
    "Individual Differences: Questionnaire Associations Replicating Across Studies 4 and 5",
    make_st9_id(st9_id_df),
    paste0(
      "Pearson correlations between person-level task features and questionnaire scores, ",
      "restricted to the five task \u00d7 questionnaire pairs showing directional replication ",
      "across Studies 4 and 5 (same sign and nominally significant in Study 5). ",
      "Task features: \u2018Mean detection confidence\u2019 = mean trial-level confidence rating; ",
      "\u2018Arousal slope: hit trials\u2019 = within-person regression slope of arousal on ",
      "breathing change magnitude, hit trials only; ",
      "\u2018Confidence\u2013accuracy correlation\u2019 = within-person r(confidence, accuracy). ",
      "All correlations are uncorrected for multiple comparisons and should be treated as exploratory. ",
      "MAIA = Multidimensional Assessment of Interoceptive Awareness (total and subscale scores). ",
      "95% CIs based on Fisher r-to-z transformation."
    )
  )

out_path <- file.path(OUTPUT_DIR, "BCAT_Supplementary_Tables.docx")
print(doc, target = out_path)
message(sprintf("Saved: %s", basename(out_path)))