# Regime Comparison Figure (Supplementary)
# Forest plot: Hit slopes | Miss slopes
# across three exclusion levels: R1 (≥3 trials), R2 (≥6 trials),
# R2b (matched magnitude), shown for all four studies
#
# Cohen's d with 95% CI from one-sample t-tests on per-person slopes
# (b_change_M1, primary linear model)
# Studies: 1A, 2, 4, 5

# Set Up ---------

# ============================================================
# 1. LOAD DATA
# ============================================================
regime_raw  <- readr::read_csv(file.path(RESULTS_DIR, "tce_sensitivity_regime.csv"))
matched_raw <- readr::read_csv(file.path(RESULTS_DIR, "tce_sensitivity_matched.csv"))

# Keep b_change_M1 only (primary linear slope)
regime <- regime_raw |>
  dplyr::filter(coef == "b_change_M1") |>
  dplyr::mutate(
    regime_label = dplyr::case_when(
      min_trials == 3 ~ "R1: ≥3 trials (primary)",
      min_trials == 6 ~ "R2: ≥6 trials",
      TRUE            ~ as.character(min_trials)
    )
  )

matched <- matched_raw |>
  dplyr::filter(coef == "b_change_M1") |>
  dplyr::mutate(
    regime_label = "R2b: matched magnitude",
    min_trials   = NA_real_
  )

# Combine and set factor order (most restrictive at top in forest plot)
plot_data <- dplyr::bind_rows(regime, matched) |>
  dplyr::mutate(
    study_label = dplyr::recode(study,
                                Study1A = "Study 1A",
                                Study2  = "Study 2",
                                Study4  = "Study 4",
                                Study5  = "Study 5"
    ),
    study_label  = factor(study_label,
                          levels = c("Study 1A", "Study 2", "Study 4", "Study 5")),
    regime_label = factor(regime_label,
                          levels = c("R2b: matched magnitude",
                                     "R2: ≥6 trials",
                                     "R1: ≥3 trials (primary)")),
    condition    = factor(condition, levels = c("hits", "misses"))
  )

# Compute d CIs from the t-distribution:
#   SE(d) = sqrt((1/n) + (d^2 / (2*n)))   [one-sample non-central t approximation]
#   CI = d ± t_{df, .975} * SE(d)
plot_data <- plot_data |>
  dplyr::mutate(
    se_d  = sqrt((1 / n_pp) + (cohens_dz^2 / (2 * n_pp))),
    t_crit = qt(.975, df = n_pp - 1),
    d_lo  = cohens_dz - t_crit * se_d,
    d_hi  = cohens_dz + t_crit * se_d
  )

# ============================================================
# 2. COLOUR PALETTE AND THEME
# ============================================================
COL_HIT  <- "#DD8452"   # orange -- hits / detected
COL_MISS <- "#4C72B0"   # blue   -- misses / undetected

# Study shapes
study_shapes <- c("Study 1A" = 15, "Study 2" = 17, "Study 4" = 16, "Study 5" = 18)
study_sizes  <- c("Study 1A" = 3,  "Study 2" = 3,  "Study 4" = 3,  "Study 5" = 3.5)

theme_forest <- function() {
  ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      strip.background    = ggplot2::element_blank(),
      strip.text          = ggplot2::element_text(face = "bold", size = 10),
      axis.title.y        = ggplot2::element_blank(),
      axis.text.y         = ggplot2::element_text(size = 9, colour = "grey30"),
      axis.text.x         = ggplot2::element_text(size = 9),
      axis.title.x        = ggplot2::element_text(size = 10),
      panel.grid.major.x  = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      legend.position     = "bottom",
      legend.title        = ggplot2::element_text(size = 9),
      legend.text         = ggplot2::element_text(size = 9),
      legend.key.size     = ggplot2::unit(0.4, "cm"),
      plot.title          = ggplot2::element_text(size = 11, face = "bold"),
      panel.spacing.y     = ggplot2::unit(0.6, "lines")
    )
}

# ============================================================
# 3. BUILD PANELS
# ============================================================
make_panel <- function(data, fill_col, title_str, x_limits, show_legend = FALSE) {
  
  p <- data |>
    ggplot2::ggplot(ggplot2::aes(
      x     = cohens_dz,
      xmin  = d_lo,
      xmax  = d_hi,
      y     = regime_label,
      shape = study_label
    )) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype   = "dashed",
      colour     = "grey40",
      linewidth  = 0.5
    ) +
    ggplot2::geom_errorbarh(
      colour = fill_col,
      height = 0.25,
      linewidth = 0.65,
      alpha  = 0.8
    ) +
    ggplot2::geom_point(
      colour = fill_col,
      fill   = fill_col,
      size   = 3,
      alpha  = 0.95
    ) +
    ggplot2::scale_shape_manual(
      name   = "Study",
      values = study_shapes
    ) +
    ggplot2::facet_wrap(~ study_label, ncol = 2, scales = "free_y") +
    ggplot2::labs(
      x     = "Cohen's d",
      title = title_str
    ) +
    ggplot2::coord_cartesian(xlim = x_limits) +
    theme_forest()
  
  if (!show_legend) p <- p + ggplot2::theme(legend.position = "none")
  p
}

# Shared x-axis limits: wide enough to show the full hit CI range,
# so that miss slopes (near zero on this scale) are visually comparable.
# The gap between panels IS the interaction story.
SHARED_XLIM <- c(-1.6, 0.7)

p_hits <- make_panel(
  data        = dplyr::filter(plot_data, condition == "hits"),
  fill_col    = COL_HIT,
  title_str   = "Detected trials: Change \u2192 Arousal slope",
  x_limits    = SHARED_XLIM,
  show_legend = FALSE
)

p_misses <- make_panel(
  data        = dplyr::filter(plot_data, condition == "misses"),
  fill_col    = COL_MISS,
  title_str   = "Missed trials: Change \u2192 Arousal slope (same scale)",
  x_limits    = SHARED_XLIM,
  show_legend = TRUE
)

# ============================================================
# 4. COMBINE WITH ANNOTATION
# ============================================================
p_combined <- (p_hits / p_misses) +
  patchwork::plot_annotation(
    title    = "Figure S[X]. Sensitivity to exclusion regime",
    subtitle = paste0(
      "Cohen's d \u00b1 95% CI for Change \u2192 Arousal slope across three exclusion criteria.\n",
      "R1 = primary (≥3 trials/cell); R2 = stricter (≥6 trials/cell); ",
      "R2b = matched magnitude restriction. Both panels use identical x-axis scale.\n",
      "Hit slopes (top, orange) are large and stable. ",
      "Miss slopes (bottom, blue) remain near zero on the same scale, ",
      "ruling out positive arousal transfer on undetected trials."
    ),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey30",
                                            lineheight = 1.4)
    )
  )

# ============================================================
# 5. SAVE
# ============================================================
ggplot2::ggsave(
  file.path(FIG_DIR, "figure_S_regime_comparison.pdf"),
  p_combined,
  width  = 8,
  height = 9,
  device = "pdf"
)

ggplot2::ggsave(
  file.path(FIG_DIR, "figure_S_regime_comparison.png"),
  p_combined,
  width  = 8,
  height = 9,
  dpi    = 300
)

cat("Saved: figure_S_regime_comparison.pdf/.png\n")
cat("Location:", FIG_DIR, "\n")