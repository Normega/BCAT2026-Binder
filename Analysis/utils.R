# ============================================================
# utils.R
# Interoception 2025 — Shared Utility Functions
#
# Source this file at the top of every analysis script:
#   source(file.path(base_dir, "Analysis/utils.R"))
#
# Contents:
#   U0.  Standard paths (set once here; used everywhere)
#   U1.  standardise_study_data()   — structural annotation per study
#   U2.  standardise_maia()         — within-study MAIA z-scoring
#   U3.  compute_awareness_index()  — within-person cor(Confidence, Accuracy)
#   U6.  load_all_data()            — load and validate all clean CSVs
#   U7.  report_lrt()               — compact LRT console output
#   U8.  extract_h4b()              — extract Change × Accuracy coefficient
#   U9.  compute_retest_icc()       — H6 test-retest ICC
#   U10. fmt_cor()                  — format a cor.test result for inline reporting
#   U11. fmt_lmer()                 — format a single lmer coefficient for inline reporting
# ============================================================


# ============================================================
# U1. standardise_study_data()
#
# Column renames now happen in individual clean scripts, so this
# function's job is narrowed to structural annotation: coercing
# id to integer, flagging study-specific quirks, and attaching
# a study_id column so downstream code always knows which study
# a row came from.
#
# Called once per study immediately after loading the CSV.
# Returns the input data frame with study_id attached and any
# structural coercions applied.
# ============================================================
standardise_study_data <- function(df, study) {
  study <- as.integer(study)

  df <- df |> dplyr::mutate(study_id = study)

  # Coerce id to integer for all studies.
  # Studies 1 & 2 use Prolific string IDs — hash to integer via factor().
  # Studies 3-5 use numeric IDs — cast directly to preserve original values.
  if ("id" %in% names(df)) {
    if (!is.integer(df$id)) {
      if (is.numeric(df$id)) {
        df <- df |> dplyr::mutate(id = as.integer(id))
      } else {
        df <- df |> dplyr::mutate(id = as.integer(factor(id)))
      }
    }
  }


  df
}


# ============================================================
# U2. standardise_maia()
#
# Z-scores MAIA total and subscales; adds _z columns.
# ============================================================
standardise_maia <- function(df, study_label = "") {
  maia_cols <- c(
    "MAIA_total",
    "MAIA_Noticing", "MAIA_NotDistracting", "MAIA_NotWorrying",
    "MAIA_AttentionReg", "MAIA_EmoAware", "MAIA_SelfReg",
    "MAIA_BodyListen", "MAIA_Trusting"
  )

  # Only z-score columns that actually exist in this study's data
  present_cols <- intersect(maia_cols, names(df))
  missing_cols <- setdiff(maia_cols, names(df))

  if (length(missing_cols) > 0 && study_label != "") {
    message(sprintf("[%s] standardise_maia: columns not present (skipped): %s",
                    study_label, paste(missing_cols, collapse = ", ")))
  }

  for (col in present_cols) {
    z_col <- paste0(col, "_z")
    df[[z_col]] <- as.numeric(scale(df[[col]]))
  }

  df
}


# ============================================================
# U3. compute_awareness_index()
#
# Within-person r(Confidence, Accuracy). Returns one row per participant.
# ============================================================
compute_awareness_index <- function(long_data, group_cols = NULL) {
  grp_vars <- c("id", group_cols)

  long_data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grp_vars))) |>
    dplyr::summarise(
      Awareness = suppressWarnings(
        cor(as.numeric(Confidence), as.numeric(Accuracy),
            use = "complete.obs")
      ),
      n_trials = dplyr::n(),
      .groups = "drop"
    ) |>
    # Flag participants with too few trials for reliable Awareness estimate
    dplyr::mutate(
      Awareness = dplyr::if_else(n_trials < 5, NA_real_, Awareness)
    )
}


# ============================================================
# U6. load_all_data()
#
# Loads all clean CSVs. Returns named list.
# ============================================================
load_all_data <- function(data_dir = DATA_DIR) {

  message("Loading clean data from: ", data_dir)

  out <- list(
    s1_long    = readr::read_csv(file.path(data_dir, "study1_long.csv")),
    s1_summary = readr::read_csv(file.path(data_dir, "study1_summary.csv")),
    s2_long    = readr::read_csv(file.path(data_dir, "study2_long.csv")),
    s2_summary = readr::read_csv(file.path(data_dir, "study2_summary.csv")),
    s3_long    = readr::read_csv(file.path(data_dir, "study3_long.csv")),
    s3_summary = readr::read_csv(file.path(data_dir, "study3_summary.csv")),
    s4_long    = readr::read_csv(file.path(data_dir, "study4_long.csv")),
    s4_test    = readr::read_csv(file.path(data_dir, "study4_test.csv")),
    s4_summary = readr::read_csv(file.path(data_dir, "study4_summary.csv")),
    s5_long    = readr::read_csv(file.path(data_dir, "study5_long.csv")),
    s5_test    = readr::read_csv(file.path(data_dir, "study5_test.csv")),
    s5_summary = readr::read_csv(file.path(data_dir, "study5_summary.csv"))
  )

  # Study 5 belt QC
  qc_path <- file.path(data_dir, "qcFile.xlsx")
  if (file.exists(qc_path)) {
    out$s5_qcFull    <- readxl::read_excel(qc_path, sheet = "FullResults")
    out$s5_qcSummary <- readxl::read_excel(qc_path, sheet = "ResultsSummary")
    message("  Belt QC loaded: FullResults (", nrow(out$s5_qcFull),
            " rows), ResultsSummary (", nrow(out$s5_qcSummary), " rows)")
  } else {
    warning("qcFile.xlsx not found at: ", qc_path)
  }

  # HBD data
  hbd_path <- file.path(data_dir, "heartbeat_detection_qc.csv")
  if (file.exists(hbd_path)) {
    out$s5_hbd <- readr::read_csv(hbd_path) |>
      dplyr::mutate(
        id          = as.integer(id),
        hbd_quality = dplyr::case_when(
          n_intervals_analyzed == 3 ~ "complete",
          n_intervals_analyzed == 0 ~ "unavailable",
          TRUE                      ~ "partial"
        ),
        # Flag three participants with implausible HR in one interval
        hbd_implausible = id %in% c(12043L, 15403L, 6313L)
      )
    message("  HBD data loaded: ", nrow(out$s5_hbd), " participants")
  } else {
    warning("heartbeat_detection_qc.csv not found at: ", hbd_path)
  }

  hbd_int_path <- file.path(data_dir, "heartbeat_detection_results.xlsx")
  if (file.exists(hbd_int_path)) {
    out$s5_hbd_intervals <- readxl::read_excel(hbd_int_path,
                                                sheet = "IntervalResults") |>
      dplyr::mutate(id = as.integer(id))
    message("  HBD intervals loaded: ", nrow(out$s5_hbd_intervals), " rows")
  } else {
    warning("heartbeat_detection_results.xlsx not found at: ", hbd_int_path)
  }

  # ── Validation checks ───────────────────────────────────────
  expected_ns <- list(
    s1_long = 4875L, s1_summary = 181L,
    s2_long = 4155L, s2_summary = 166L,
    s3_long = 12768L, s3_summary = 103L,
    s4_long = 5240L,  s4_summary = 131L,
    s5_long = 14520L, s5_summary = 206L
  )

  all_ok <- TRUE
  for (nm in names(expected_ns)) {
    if (nm %in% names(out)) {
      actual   <- nrow(out[[nm]])
      expected <- expected_ns[[nm]]
      if (actual != expected) {
        warning(sprintf(
          "Row count mismatch for %s: expected %d, got %d. Re-run clean script?",
          nm, expected, actual
        ))
        all_ok <- FALSE
      } else {
        message(sprintf("  %-14s %d rows OK", nm, actual))
      }
    }
  }
  if (all_ok) message("All row counts validated.")

  out
}


# ============================================================
# U7. report_lrt()
#
# Prints a compact LRT table with R² for each merMod to console.
# ============================================================
report_lrt <- function(lrt_object, study_label = "") {
  cat(sprintf("\n=== LRT: %s ===\n", study_label))
  print(lrt_object)

  # R² for each model in the comparison sequence
  mods <- attr(lrt_object, "models")
  if (!is.null(mods)) {
    cat("R² (marginal / conditional):\n")
    for (m in mods) {
      if (inherits(m, "merMod")) {
        r2 <- tryCatch(
          MuMIn::r.squaredGLMM(m),
          error = function(e) c(NA_real_, NA_real_)
        )
        cat(sprintf("  R2m = %.3f  R2c = %.3f\n", r2[1], r2[2]))
      }
    }
  }
  invisible(lrt_object)
}


# ============================================================
# U8. extract_h4b()
#
# Extracts the Change × Accuracy interaction coefficient from a
# fitted lmer/lmerMod object. Returns a named numeric vector
# c(b, se, t, df, p, partial_r) or NAs if the term is not found.
#
# df is taken from the lmerTest Satterthwaite column when available;
# falls back to nobs - nfixef so partial_r is always exact.
#
# Used to populate the cross-study H4B summary table.
# ============================================================
extract_h4b <- function(model, study_label = "") {
  na_out <- c(b = NA_real_, se = NA_real_, t = NA_real_,
              df = NA_real_, p = NA_real_, partial_r = NA_real_)

  if (!inherits(model, "merMod")) {
    message(sprintf("[%s] extract_h4b: not a merMod object", study_label))
    return(na_out)
  }

  cf <- summary(model)$coefficients
  # Match the interaction term regardless of order
  int_row <- grep("Change.*Accuracy|Accuracy.*Change",
                  rownames(cf), value = TRUE)

  if (length(int_row) == 0) {
    message(sprintf("[%s] extract_h4b: Change × Accuracy term not found",
                    study_label))
    return(na_out)
  }

  row <- cf[int_row[1], , drop = FALSE]
  b   <- row[1, "Estimate"]
  se  <- row[1, "Std. Error"]
  t   <- row[1, "t value"]
  # lmerTest adds a "df" column (Satterthwaite); lme4::lmer does not.
  # Fall back to a residual df approximation so partial_r is always exact.
  df  <- if ("df" %in% colnames(cf)) row[1, "df"] else
           nobs(model) - length(lme4::fixef(model)) -
           sum(sapply(lme4::ranef(model), nrow))
  p   <- row[1, "Pr(>|t|)"]
  pr  <- t / sqrt(t^2 + df)

  c(b = b, se = se, t = t, df = df, p = p, partial_r = pr)
}


# ============================================================
# U9. compute_retest_icc()
#
# Two-way agreement ICC for test-retest reliability.
# ============================================================
compute_retest_icc <- function(summary_data, ses1_col, ses2_col,
                               label = "") {
  d <- summary_data |>
    dplyr::select(id,
                  ses1 = dplyr::all_of(ses1_col),
                  ses2 = dplyr::all_of(ses2_col)) |>
    tidyr::drop_na()

  icc_result <- irr::icc(
    cbind(d$ses1, d$ses2),
    model = "twoway", type = "agreement", unit = "single"
  )

  list(label = label, n = nrow(d), icc = icc_result)
}


# ============================================================
# U10. fmt_cor()
#
# Format cor.test as inline string: "r(163) = .42, p = .003"
# ============================================================
fmt_cor <- function(ct, digits = 2) {
  r   <- round(ct$estimate,  digits)
  p   <- ct$p.value
  df  <- ct$parameter
  p_s <- dplyr::if_else(p < .001, "< .001", paste0("= ", round(p, 3)))
  sprintf("r(%d) = %s, p %s",
          as.integer(df),
          formatC(r, format = "f", digits = digits),
          p_s)
}


# ============================================================
# U12. make_threshold_long()
#
# Pivot threshold columns to long format (one row per participant × condition).
# Returns NULL for Studies 2 and 3.
# ============================================================
make_threshold_long <- function(summary_df, study,
                                ses_filter = NULL) {
  if (study == 1) {
    summary_df |>
      dplyr::select(id, thresh_TaskA, thresh_TaskB) |>
      tidyr::pivot_longer(
        cols      = c(thresh_TaskA, thresh_TaskB),
        names_to  = "Group",
        values_to = "Threshold"
      ) |>
      dplyr::mutate(
        Group     = stringr::str_remove(Group, "thresh_"),
        Salience  = NA_character_,
        Direction = NA_character_
      )
  } else if (study %in% c(2, 3)) {
    NULL   # no staircase threshold in Study 2/3 summary
  } else if (study == 4) {
    summary_df |>
      tidyr::pivot_longer(
        cols          = dplyr::matches("^thresh_"),
        names_to      = c("Salience", "Direction"),
        names_pattern = "thresh_([^_]+)_(.*)",
        values_to     = "Threshold"
      ) |>
      dplyr::mutate(
        Salience  = factor(Salience,  levels = c("Low", "High")),
        Direction = factor(Direction, levels = c("Faster", "Slower"))
      )
  } else if (study == 5) {
    long_df <- summary_df |>
      tidyr::pivot_longer(
        cols          = dplyr::matches("^thresh_ses[12]_"),
        names_to      = c("ses", "Salience", "Direction"),
        names_pattern = "thresh_(ses[12])_([^_]+)_(.*)",
        values_to     = "Threshold"
      ) |>
      dplyr::mutate(
        Salience  = factor(Salience,  levels = c("Low", "High")),
        Direction = factor(Direction, levels = c("Faster", "Slower")),
        ses       = factor(ses)
      )
    if (!is.null(ses_filter))
      long_df <- dplyr::filter(long_df, ses %in% ses_filter)
    long_df
  }
}


# ============================================================
# U11. fmt_lmer()
#
# Format lmerTest fixed effect as inline string.
# ============================================================
fmt_lmer <- function(model, term, digits = 2) {
  cf <- summary(model)$coefficients
  if (!term %in% rownames(cf)) {
    warning(sprintf("fmt_lmer: term '%s' not found in model", term))
    return(NA_character_)
  }
  row <- cf[term, , drop = FALSE]
  b   <- round(row[1, "Estimate"],   digits)
  se  <- round(row[1, "Std. Error"], digits)
  tv  <- round(row[1, "t value"],    2)

  # Pr(>|t|) only present when lmerTest::lmer() was used
  p_s <- tryCatch({
    p <- row[1, "Pr(>|t|)"]
    dplyr::if_else(p < .001, "< .001", paste0("= ", round(p, 3)))
  }, error = function(e) "= NA (refit with lmerTest::lmer)")

  df_val <- tryCatch(round(row[1, "df"], 0), error = function(e) "")

  sprintf("b = %s, SE = %s, t(%s) = %s, p %s",
          formatC(b,  format = "f", digits = digits),
          formatC(se, format = "f", digits = digits),
          df_val, formatC(tv, format = "f", digits = 2),
          p_s)
}


# ── Effect size helpers (used across all analysis files) ─────
# ── Effect size helpers ───────────────────────────────────────

# Partial r from LMM/GLMM fixed effect: t / sqrt(t^2 + df_resid)
# When df unknown, recover from two-sided p-value via uniroot.
partial_r_from_t <- function(b, se, p = NULL, df = NULL) {
  if (is.na(b) || is.na(se)) return(NA_real_)
  t_val <- b / se
  if (!is.null(df) && !is.na(df))
    return(t_val / sqrt(t_val^2 + df))
  if (!is.null(p) && !is.na(p) && p > 0 && p < 1) {
    df_est <- tryCatch(
      stats::uniroot(
        function(df) 2 * stats::pt(-abs(t_val), df) - p,
        lower = 1, upper = 1e6)$root,
      error = function(e) NA_real_)
    if (!is.na(df_est))
      return(t_val / sqrt(t_val^2 + df_est))
  }
  NA_real_
}

# Odds ratio + 95% CI from GLMM log-odds coefficient
odds_ratio <- function(b, se) {
  tibble::tibble(
    OR       = exp(b),
    OR_lower = exp(b - 1.96 * se),
    OR_upper = exp(b + 1.96 * se)
  )
}

# Add OR columns to a data frame (expects b_ and se_ prefixes)
add_odds_ratios <- function(df, terms = c("Change", "Change2")) {
  for (term in terms) {
    b_col  <- paste0("b_",  term)
    se_col <- paste0("se_", term)
    if (!all(c(b_col, se_col) %in% names(df))) next
    df[[paste0("OR_", term)]]           <- exp(df[[b_col]])
    df[[paste0("OR_", term, "_lower")]] <- exp(df[[b_col]] - 1.96 * df[[se_col]])
    df[[paste0("OR_", term, "_upper")]] <- exp(df[[b_col]] + 1.96 * df[[se_col]])
  }
  df
}

# Add partial_r columns to a data frame (expects b_ se_ p_ prefixes)
add_partial_r <- function(df, hypotheses) {
  for (hyp in hypotheses) {
    b_col <- paste0(hyp, "_b")
    s_col <- paste0(hyp, "_SE")
    p_col <- paste0(hyp, "_p")
    if (!all(c(b_col, s_col) %in% names(df))) next
    df[[paste0(hyp, "_partial_r")]] <- mapply(
      partial_r_from_t,
      b  = df[[b_col]],
      se = df[[s_col]],
      p  = if (p_col %in% names(df)) df[[p_col]] else rep(NA_real_, nrow(df))
    )
  }
  df
}

# Extract Change and Change² betas + LRT from fit_detection_models() output
extract_change2_results <- function(det_obj, study_label, n_participants) {
  cf  <- summary(det_obj$quadratic)$coefficients
  lrt <- det_obj$lrt_main
  tibble::tibble(
    study        = study_label,
    n            = n_participants,
    b_Change     = cf["Change",  "Estimate"],
    se_Change    = cf["Change",  "Std. Error"],
    z_Change     = cf["Change",  "z value"],
    p_Change     = cf["Change",  "Pr(>|z|)"],
    b_Change2    = cf["Change2", "Estimate"],
    se_Change2   = cf["Change2", "Std. Error"],
    z_Change2    = cf["Change2", "z value"],
    p_Change2    = cf["Change2", "Pr(>|z|)"],
    lrt_p_linear = tryCatch(lrt[["Pr(>Chisq)"]][2], error = function(e) NA_real_),
    lrt_p_quad   = tryCatch(lrt[["Pr(>Chisq)"]][3], error = function(e) NA_real_)
  )
}

# Group-level 3AFC d' from test-block Pc, via numerical inversion.
compute_test_dprime_3afc <- function(test_data,
                                      study_label,
                                      group_filter = NULL,
                                      group_col    = "Group",
                                      acc_col      = "Accuracy",
                                      salience_col = "Salience",
                                      id_col       = "id") {
  if (!is.null(group_filter))
    test_data <- dplyr::filter(test_data, .data[[group_col]] == group_filter)

  pc_from_dprime <- function(d) {
    stats::integrate(
      function(x) stats::dnorm(x - d) * stats::pnorm(x)^2,
      lower = -10, upper = 10)$value
  }
  dprime_from_pc <- function(pc) {
    if (pc <= 1/3) return(0)
    upper <- 1
    while (pc_from_dprime(upper) < pc) {
      upper <- upper + 1
      if (upper > 20) return(20)
    }
    stats::uniroot(function(d) pc_from_dprime(d) - pc,
                   lower = 0, upper = upper)$root
  }

  purrr::map_dfr(c("High", "Low"), function(sal) {
    sub <- dplyr::filter(test_data, .data[[salience_col]] == sal)
    pc  <- mean(sub[[acc_col]], na.rm = TRUE)
    tibble::tibble(
      study          = study_label,
      group          = group_filter %||% "All",
      salience       = sal,
      n_trials       = nrow(sub),
      n_participants = dplyr::n_distinct(sub[[id_col]]),
      Pc             = round(pc, 3),
      dprime_3afc    = round(dprime_from_pc(pc), 3)
    )
  })
}

`%||%` <- function(a, b) if (!is.null(a)) a else b


# ── BH correction ────────────────────────────────────────────
# Applies Benjamini-Hochberg correction to all 8 pre-registered hypothesis
# p-values for a single study, treated as one correction family.
# Correct Studies 4 and 5 separately by calling this function twice.
#
# Expected hypothesis names: H1, H2, H3A, H3B, H4A, H5, H4B, H4C
#
# pvals must be a named numeric vector with names exactly matching the
# hypothesis labels above. Names like "H1.p" indicate a single-bracket
# indexing bug in the caller — use [["p"]] not ["p"] when extracting from
# named vectors returned by .extract_coef() etc.
apply_bh_correction <- function(pvals, study_label = "") {

  all_hypotheses <- c("H1", "H2", "H3A", "H3B", "H4A", "H5", "H4B", "H4C")

  # Defensive check: warn if ".p"-suffixed names are found, indicating a
  # single-bracket indexing bug in the caller.
  dotted <- paste0(all_hypotheses, ".p")
  found_dotted <- dotted[dotted %in% names(pvals)]
  if (length(found_dotted) > 0) {
    warning(sprintf(
      "[%s] apply_bh_correction: found names %s in pvals. ",
      study_label, paste(found_dotted, collapse = ", "),
      "Use [['p']] not ['p'] when extracting from named vectors. ",
      "Affected hypotheses will be treated as NA."
    ))
  }

  result <- tibble::tibble(
    hypothesis = all_hypotheses,
    p_raw      = pvals[all_hypotheses]
  ) |>
    dplyr::filter(!is.na(p_raw)) |>
    dplyr::mutate(
      p_BH   = p.adjust(p_raw, method = "BH"),
      sig_BH = p_BH < .05
    )

  cat(sprintf("\n[%s] BH correction results (all 8 hypotheses):\n", study_label))
  print(as.data.frame(result), digits = 4, row.names = FALSE)

  result
}