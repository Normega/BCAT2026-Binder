# ============================================================
# meta_analysis.R
# Random-effects meta-analytic synthesis across studies.
# Sourced by MainAnalysis.R after per-study models are fitted.
#
# Functions:
#   run_h4b_meta()      — H4B: Change x Accuracy on Arousal (k=5)
#   run_h3_meta()       — H3A/H3B: MAIA dissociation
#   run_h1_meta()       — H1: direction threshold (defined, not yet called)
#   run_h2_meta()       — H2: salience threshold (defined, not yet called)
#   write_meta_results() — save CSVs
#   plot_meta_figures()  — forest plots
#
# Outputs (CSVs):
#   meta_h4b_study_estimates.csv   — per-study H4B partial r
#   meta_h4b_pooled.csv            — pooled H4B random-effects estimate
#   meta_h3_maia_dissociation.csv  — MAIA confidence vs threshold meta
#   meta_h3_confidence_studies.csv — per-study MAIA-confidence estimates
#   meta_h3_threshold_studies.csv  — per-study MAIA-threshold estimates
#   meta_h1_direction_studies.csv  — per-study direction threshold estimates
#   meta_h2_salience_studies.csv   — per-study salience threshold estimates
#
# Outputs (figures):
#   MetaAnalysis/meta_h4b_forest.pdf          — H4B forest plot
#   MetaAnalysis/meta_h3_maia_dissociation.pdf — MAIA dissociation forest plot
# ============================================================

meta_fig_dir <- file.path(FIG_DIR, "MetaAnalysis")
dir.create(meta_fig_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================
# Internal helper: extract_study_h4b()
#
# Pulls Change × Accuracy coefficient + SE from a fitted lmer
# H4B model and converts to partial r for cross-study comparability.
# t_to_r: r = t / sqrt(t^2 + df_resid)
# SE of partial r: sqrt((1-r^2)^2 / (df_resid-1)) — uses residual df not n
# ============================================================
.extract_study_h4b <- function(model, n, study_label) {
  cf   <- summary(model)$coefficients
  term <- grep("Change.*Accuracy|Accuracy.*Change",
               rownames(cf), value = TRUE)
  if (length(term) == 0) {
    warning(sprintf("[%s] Change×Accuracy term not found", study_label))
    return(NULL)
  }
  row  <- cf[term[1], , drop = FALSE]
  b    <- row[1, "Estimate"]
  se   <- row[1, "Std. Error"]
  tval <- row[1, "t value"]
  p    <- tryCatch(row[1, "Pr(>|t|)"], error = function(e) NA_real_)
  df_r <- tryCatch(row[1, "df"],       error = function(e) n - 3)

  r_partial <- tval / sqrt(tval^2 + df_r)
  # SE of partial r: use residual df from the multilevel model, not
  # participant n. The model df already accounts for the nested structure;
  # using n would massively inflate SE relative to the model's actual precision.
  se_r      <- sqrt((1 - r_partial^2)^2 / (df_r - 1))

  tibble::tibble(
    study     = study_label,
    n         = n,
    b_raw     = b,
    se_raw    = se,
    t_val     = tval,
    df_resid  = df_r,
    r_partial = r_partial,
    se_r      = se_r,
    p_value   = p
  )
}


# ============================================================
# MA1. run_h4b_meta()
#
# Random-effects meta-analysis of Change × Accuracy on Arousal.
# `study_results`: named list; each element must have $H4B_final.
# `study_ns`:      named list of sample sizes matching study_results.
# ============================================================
run_h4b_meta <- function(study_results, study_ns) {
  cat("\n\n========================================\n")
  cat("MA1: H4B — Change × Accuracy on Arousal\n")
  cat("Random-effects meta-analysis (partial r)\n")
  cat("========================================\n")

  h4b_rows <- purrr::map_dfr(names(study_results), function(nm) {
    res <- study_results[[nm]]
    n   <- study_ns[[nm]]
    if (is.null(res$H4B_final)) {
      message(sprintf("  [%s] H4B_final not found — skipping", nm))
      return(NULL)
    }
    .extract_study_h4b(res$H4B_final, n = n, study_label = nm)
  })

  if (nrow(h4b_rows) < 2) {
    warning("Fewer than 2 H4B results — meta-analysis skipped")
    return(NULL)
  }

  cat(sprintf("\nStudies included: %d\n", nrow(h4b_rows)))
  print(
    h4b_rows |>
      dplyr::select(study, n, r_partial, se_r, p_value) |>
      dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, 4))) |>
      as.data.frame(),
    row.names = FALSE
  )

  rma_out <- metafor::rma(
    yi     = r_partial,
    sei    = se_r,
    data   = h4b_rows,
    method = "REML",
    slab   = h4b_rows$study
  )

  cat("\nPooled random-effects estimate:\n")
  print(summary(rma_out))
  cat(sprintf(
    "Heterogeneity: τ²=%.4f  I²=%.1f%%  Q(%d)=%.2f  p=%.4f\n",
    rma_out$tau2, rma_out$I2,
    rma_out$k - 1, rma_out$QE, rma_out$QEp
  ))

  list(data = h4b_rows, rma = rma_out)
}


# ============================================================
# MA2. run_h3_meta()
#
# Pools two MAIA correlations separately via Fisher z:
#   r(Confidence, MAIA) — H3A, expect positive
#   r(Threshold,  MAIA) — H3B, expect null (BF₀₁ evidence)
#
# `maia_table`: tibble with columns study, n_Confidence,
#   r_Confidence_MAIA, n_Threshold, r_Threshold_MAIA.
#   Built in unified_analysis.R from individual fit_maia_tests() results.
# ============================================================
run_h3_meta <- function(maia_table) {
  cat("\n\n========================================\n")
  cat("MA2: H3 — MAIA Dissociation\n")
  cat("Pooled correlations via Fisher z\n")
  cat("========================================\n")

  results <- list()

  for (outcome in c("Confidence", "Threshold")) {
    r_col <- paste0("r_", outcome, "_MAIA")
    n_col <- paste0("n_", outcome)

    if (!all(c(r_col, n_col) %in% names(maia_table))) {
      message(sprintf("  %s columns not found — skipping", outcome))
      next
    }

    d <- maia_table |>
      dplyr::filter(!is.na(.data[[r_col]]),
                    !is.na(.data[[n_col]])) |>
      dplyr::rename(r = dplyr::all_of(r_col),
                    n = dplyr::all_of(n_col))

    if (nrow(d) < 2) next

    d_esc <- metafor::escalc(
      measure = "COR",
      ri   = d$r,
      ni   = d$n,
      data = d,
      slab = d$study
    )

    rma_out <- metafor::rma(
      yi = yi, vi = vi,
      data   = d_esc,
      method = "REML"
    )

    r_pooled <- tanh(rma_out$beta[[1]])
    r_lower  <- tanh(rma_out$ci.lb)
    r_upper  <- tanh(rma_out$ci.ub)

    cat(sprintf(
      "\nH3 %s ~ MAIA: pooled r = %.3f [%.3f, %.3f]  p = %.4f",
      outcome, r_pooled, r_lower, r_upper, rma_out$pval
    ))
    cat(sprintf(
      "  I²=%.1f%%  τ²=%.4f\n",
      rma_out$I2, rma_out$tau2
    ))

    results[[outcome]] <- list(
      data     = d_esc,
      rma      = rma_out,
      r_pooled = r_pooled,
      r_lower  = r_lower,
      r_upper  = r_upper
    )
  }

  results$summary <- tibble::tibble(
    Contrast   = c("Confidence ~ MAIA (H3A)",
                   "Threshold  ~ MAIA (H3B)"),
    k_studies  = c(
      if (!is.null(results$Confidence)) results$Confidence$rma$k else NA_integer_,
      if (!is.null(results$Threshold))  results$Threshold$rma$k  else NA_integer_
    ),
    r_pooled   = c(
      if (!is.null(results$Confidence)) results$Confidence$r_pooled else NA_real_,
      if (!is.null(results$Threshold))  results$Threshold$r_pooled  else NA_real_
    ),
    r_lower    = c(
      if (!is.null(results$Confidence)) results$Confidence$r_lower  else NA_real_,
      if (!is.null(results$Threshold))  results$Threshold$r_lower   else NA_real_
    ),
    r_upper    = c(
      if (!is.null(results$Confidence)) results$Confidence$r_upper  else NA_real_,
      if (!is.null(results$Threshold))  results$Threshold$r_upper   else NA_real_
    ),
    p_value    = c(
      if (!is.null(results$Confidence)) results$Confidence$rma$pval else NA_real_,
      if (!is.null(results$Threshold))  results$Threshold$rma$pval  else NA_real_
    ),
    # Study 3 excluded from pooled H3B only (accuracy proxy, not threshold).
    # H3A pools all 5 studies — confidence ratings available in Study 3.
    note       = c(
      "All 5 studies included",
      "Study 3 excluded (accuracy proxy, not threshold)"
    )
  )

  cat("\nMAIA dissociation — pooled summary:\n")
  print(results$summary, digits = 3)
  results
}


# ============================================================
# MA3. run_h1_meta() — Direction asymmetry on threshold
# MA4. run_h2_meta() — Salience effect on threshold
#
# Both expect an input tibble with columns: study, n, d, se_d
# (Cohen's d and its SE, extracted from threshold models in unified).
# ============================================================
run_h1_meta <- function(h1_table) {
  cat("\n\n========================================\n")
  cat("MA3: H1 — Direction Asymmetry (threshold)\n")
  cat("========================================\n")
  if (is.null(h1_table) || nrow(h1_table) < 2) {
    message("Fewer than 2 H1 estimates — skipping")
    return(NULL)
  }
  rma_out <- metafor::rma(yi = d, sei = se_d, data = h1_table,
                           method = "REML", slab = h1_table$study)
  cat("\nH1 pooled:\n"); print(summary(rma_out))
  list(data = h1_table, rma = rma_out)
}

run_h2_meta <- function(h2_table) {
  cat("\n\n========================================\n")
  cat("MA4: H2 — Salience Effect (threshold)\n")
  cat("========================================\n")
  if (is.null(h2_table) || nrow(h2_table) < 2) {
    message("Fewer than 2 H2 estimates — skipping")
    return(NULL)
  }
  rma_out <- metafor::rma(yi = d, sei = se_d, data = h2_table,
                           method = "REML", slab = h2_table$study)
  cat("\nH2 pooled:\n"); print(summary(rma_out))
  list(data = h2_table, rma = rma_out)
}


# ============================================================
# MA5. write_meta_results()
# ============================================================
write_meta_results <- function(h4b_meta, h3_meta,
                               h1_meta = NULL, h2_meta = NULL) {
  cat("\nWriting meta-analysis results...\n")

  if (!is.null(h4b_meta)) {
    readr::write_csv(
      h4b_meta$data,
      file.path(RESULTS_DIR, "meta_h4b_study_estimates.csv")
    )
    rma <- h4b_meta$rma
    tibble::tibble(
      analysis  = "H4B Change×Accuracy",
      k_studies = rma$k,
      r_pooled  = tanh(rma$beta[[1]]),
      r_lower   = tanh(rma$ci.lb),
      r_upper   = tanh(rma$ci.ub),
      p_value   = rma$pval,
      I2_pct    = rma$I2,
      tau2      = rma$tau2,
      Q_stat    = rma$QE,
      Q_p       = rma$QEp
    ) |> readr::write_csv(file.path(RESULTS_DIR, "meta_h4b_pooled.csv"))
    message("  meta_h4b_*.csv")
  }

  if (!is.null(h3_meta$summary)) {
    readr::write_csv(h3_meta$summary,
                     file.path(RESULTS_DIR, "meta_h3_maia_dissociation.csv"))
    for (out in c("Confidence", "Threshold")) {
      if (!is.null(h3_meta[[out]])) {
        readr::write_csv(
          as.data.frame(h3_meta[[out]]$data),
          file.path(RESULTS_DIR, paste0("meta_h3_", tolower(out), "_studies.csv"))
        )
      }
    }
    message("  meta_h3_*.csv")
  }

  if (!is.null(h1_meta))
    readr::write_csv(h1_meta$data,
                     file.path(RESULTS_DIR, "meta_h1_direction_studies.csv"))
  if (!is.null(h2_meta))
    readr::write_csv(h2_meta$data,
                     file.path(RESULTS_DIR, "meta_h2_salience_studies.csv"))

  message("Meta-analysis results written to: ", RESULTS_DIR)
}


# ============================================================
# MA6. plot_meta_figures()
# ============================================================
plot_meta_figures <- function(h4b_meta, h3_meta) {

    meta_theme <- list(
    theme_bcat(base_size = 12),
    ggplot2::theme(
      plot.subtitle    = ggplot2::element_text(colour = "grey40", size = 10),
      axis.line.y      = ggplot2::element_blank(),
      axis.ticks.y     = ggplot2::element_blank(),
      panel.spacing    = ggplot2::unit(2, "lines"),
      legend.position  = "none"
    )
  )

  # Internal helper: build a single forest panel from a tidy data frame.
  # `plot_df` columns: label (chr), r_val, r_lower, r_upper, is_pooled (lgl).
  # "Pooled" row must be present and named exactly "Pooled".
  .forest_gg <- function(plot_df, xlab, title, subtitle = NULL) {

    # Extract study labels preserving their arrival order, exclude Pooled
    study_labels <- as.character(plot_df$label[plot_df$label != "Pooled"])

    # Display order: Study1 at top, Pooled at bottom
    # No separator row in the data — just a geom_hline between studies and Pooled
    display_order <- c(study_labels, "Pooled")

    pd <- plot_df |> dplyr::mutate(label = as.character(label))

    # Separator drawn between Pooled (position 1 from bottom) and
    # the last study (position 2 from bottom)
    sep_y <- 1.5

    ggplot2::ggplot(
      pd,
      ggplot2::aes(x = r_val, y = label,
                   xmin = r_lower, xmax = r_upper)
    ) +
      ggplot2::scale_y_discrete(limits = rev(display_order)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                          colour = "grey60", linewidth = 0.5) +
      ggplot2::geom_hline(yintercept = sep_y,
                          linewidth = 0.5, colour = "grey30") +
      # Individual study CIs — black lines with vertical end caps
      ggplot2::geom_errorbarh(
        data   = \(d) dplyr::filter(d, !is_pooled, label != "—"),
        height = 0.2, linewidth = 0.7, colour = "black"
      ) +
      # Individual study points — black filled squares
      ggplot2::geom_point(
        data   = \(d) dplyr::filter(d, !is_pooled, label != "—"),
        shape  = 15, size = 3.5, colour = "black"
      ) +
      # Pooled CI — blue line with end caps
      ggplot2::geom_errorbarh(
        data   = \(d) dplyr::filter(d, is_pooled),
        height = 0.2, linewidth = 1.0, colour = BCAT_BLUE
      ) +
      # Pooled point — blue diamond
      ggplot2::geom_point(
        data   = \(d) dplyr::filter(d, is_pooled),
        shape  = 18, size = 5.5, colour = BCAT_BLUE
      ) +
      ggplot2::labs(title = title, subtitle = subtitle,
                    x = xlab, y = NULL) +
      meta_theme
  }

  # ── F1: H4B ──────────────────────────────────────────────────
  if (!is.null(h4b_meta)) {
    rma <- h4b_meta$rma

    f1_data <- h4b_meta$data |>
      dplyr::mutate(
        r_lower   = tanh(atanh(r_partial) - 1.96 / sqrt(df_resid - 1)),
        r_upper   = tanh(atanh(r_partial) + 1.96 / sqrt(df_resid - 1)),
        is_pooled = FALSE,
        label     = study
      ) |>
      dplyr::select(label, r_val = r_partial, r_lower, r_upper, is_pooled) |>
      dplyr::bind_rows(tibble::tibble(
        label     = "Pooled",
        r_val     = tanh(rma$beta[[1]]),
        r_lower   = tanh(rma$ci.lb),
        r_upper   = tanh(rma$ci.ub),
        is_pooled = TRUE
      ))

    f1 <- .forest_gg(
      f1_data,
      xlab     = "Partial r",
      title    = "H4B: Change x Accuracy -> Arousal",
      subtitle = sprintf(
        "Pooled r = %.3f [%.3f, %.3f]  I\u00b2 = %.0f%%  p = %.4f",
        tanh(rma$beta[[1]]), tanh(rma$ci.lb), tanh(rma$ci.ub),
        rma$I2, rma$pval
      )
    )

    ggplot2::ggsave(file.path(meta_fig_dir, "meta_h4b_forest.pdf"),
                   f1, width = 8, height = 5.5, device = "pdf")
    message("Saved: meta_h4b_forest.pdf")
  }

  # ── F2: H3 two-panel ─────────────────────────────────────────
  h3_panels <- purrr::compact(purrr::map(
    stats::setNames(c("Confidence","Threshold"), c("Confidence","Threshold")),
    function(outcome) {
      obj <- h3_meta[[outcome]]
      if (is.null(obj)) return(NULL)
      as.data.frame(obj$data) |>
        dplyr::mutate(
          r_val     = tanh(yi),
          r_lower   = tanh(yi - 1.96 * sqrt(vi)),
          r_upper   = tanh(yi + 1.96 * sqrt(vi)),
          label     = as.character(study),
          is_pooled = FALSE
        ) |>
        dplyr::select(label, r_val, r_lower, r_upper, is_pooled) |>
        dplyr::bind_rows(tibble::tibble(
          label     = "Pooled",
          r_val     = obj$r_pooled,
          r_lower   = obj$r_lower,
          r_upper   = obj$r_upper,
          is_pooled = TRUE
        ))
    }
  ))

  if (length(h3_panels) > 0) {
    panel_titles <- c(
      Confidence = "Confidence ~ MAIA (H3A: expect +)",
      Threshold  = "Threshold ~ MAIA (H3B: expect null)"
    )
    plots <- purrr::imap(h3_panels, function(pd, nm) {
      rma_obj <- h3_meta[[nm]]$rma
      .forest_gg(
        pd,
        xlab     = "Correlation (r)",
        title    = panel_titles[[nm]],
        subtitle = sprintf(
          "Pooled r = %.3f [%.3f, %.3f]  I\u00b2 = %.0f%%",
          h3_meta[[nm]]$r_pooled,
          h3_meta[[nm]]$r_lower,
          h3_meta[[nm]]$r_upper,
          rma_obj$I2
        )
      )
    })

    f2 <- patchwork::wrap_plots(plots, nrow = 1) +
      patchwork::plot_annotation(
        title = "H3: MAIA dissociation across studies",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 13)
        )
      )

    ggplot2::ggsave(file.path(meta_fig_dir, "meta_h3_maia_dissociation.pdf"),
                   f2, width = 12, height = 5.5, device = "pdf")
    message("Saved: meta_h3_maia_dissociation.pdf")
  }
}
