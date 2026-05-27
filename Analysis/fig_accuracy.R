# fig_accuracy.R
# Change magnitude predicts detection accuracy across studies
# Layout mirrors fig_replication.R: 3 columns × 2 rows
#   Row 1: Study 1A | Study 2 | Study 3
#   Row 2: Study 4  | Study 5 | Forest plot (Change -> Accuracy partial r)
#
# Studies 3, 4 and 5 split by Salience (High / Low).
# Studies 1A and 2 have no salience manipulation — single smooth.
# Forest plot: partial r of Change on Accuracy per study (GLM-derived).

# Set Up ---------

ctrl     <- lmerControl(optimizer = "bobyqa")

# ── Colour palette ──────────────────────────────────────────────────────────
# High salience = orange (same as "Detected" in fig_replication)
# Low  salience = blue   (same as "Missed"   in fig_replication)
col_high <- "#E69F00"
col_low  <- "#0072B2"

sal_cols   <- c("High" = col_high, "Low" = col_low)
sal_labels <- c("High" = "High salience", "Low" = "Low salience")

# ── Shared theme ────────────────────────────────────────────────────────────
theme_fig <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 11),
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      panel.grid.minor = element_blank(),
      legend.position  = "none"
    )
}

# ── Helper: partial r for a single named predictor from glmer ───────────────
# Fits glmer(Accuracy ~ pred + (1|id), family=binomial) and returns
# partial r with 95% CI using t/z to partial-r conversion.
# Uses approximate df = n_participants - 2.

.partial_r_glmer <- function(df, pred_name) {
  df <- dplyr::filter(df, !is.na(Accuracy), !is.na(.data[[pred_name]]))
  if (nrow(df) < 20) return(tibble::tibble(r=NA, lower=NA, upper=NA))

  m <- tryCatch(
    lme4::glmer(
      stats::as.formula(paste("Accuracy ~", pred_name, "+ (1 | id)")),
      data    = df,
      family  = binomial,
      control = lme4::glmerControl(optimizer = "bobyqa")
    ),
    error = function(e) NULL
  )
  if (is.null(m)) return(tibble::tibble(r=NA, lower=NA, upper=NA))

  n_pp  <- dplyr::n_distinct(df$id)
  df_r  <- n_pp - 2L
  b     <- lme4::fixef(m)[[pred_name]]
  s     <- sqrt(lme4::vcov.merMod(m)[pred_name, pred_name])
  z     <- b / s

  r_val  <- z / sqrt(z^2 + df_r)
  z_lo   <- (b - 1.96 * s) / s
  z_hi   <- (b + 1.96 * s) / s
  r_lo   <- z_lo / sqrt(z_lo^2 + df_r)
  r_hi   <- z_hi / sqrt(z_hi^2 + df_r)

  tibble::tibble(r = round(r_val, 3),
                 lower = round(r_lo, 3),
                 upper = round(r_hi, 3))
}

# Compute partial r for |Change| and Change² per study
.study_forest_rows <- function(df, study_label) {
  df <- df |>
    dplyr::filter(Direction %in% c("Faster", "Slower"),
                  !is.na(Accuracy), !is.na(Change)) |>
    dplyr::mutate(abs_change  = abs(Change),
                  change_sq   = Change^2)

  r_abs <- .partial_r_glmer(df, "abs_change")
  r_sq  <- .partial_r_glmer(df, "change_sq")

  dplyr::bind_rows(
    dplyr::mutate(r_abs, study = study_label, term = "|Change|"),
    dplyr::mutate(r_sq,  study = study_label, term = "Change\u00b2")
  )
}

# ── Panel builder: no salience split ────────────────────────────────────────
make_acc_panel <- function(df, title, subtitle = NULL,
                            show_y = TRUE,
                            x_breaks = c(-0.5, 0, 0.5),
                            x_labels = c("-0.5", "0", "0.5")) {
  df <- dplyr::filter(df,
                      Direction %in% c("Faster", "Slower"),
                      !is.na(Accuracy), !is.na(Change))

  m <- tryCatch(
    lme4::glmer(Accuracy ~ poly(Change, 2) + (1 | id),
                data = df, family = binomial,
                control = glmerControl(optimizer = "bobyqa")),
    error = function(e)
      lme4::glmer(Accuracy ~ Change + (1 | id),
                  data = df, family = binomial,
                  control = glmerControl(optimizer = "bobyqa"))
  )

  pred <- ggeffects::ggpredict(m, terms = "Change [all]")

  ggplot(pred, aes(x = x, y = predicted)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                fill = col_high, alpha = 0.25) +
    geom_line(colour = col_high, linewidth = 1) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = c(0.29, 1),
                       breaks = c(0.33, 0.67, 1),
                       labels = c(".33", ".67", "1")) +
    labs(x        = "Breathing rate change",
         y        = if (show_y) "P(correct detection)" else NULL,
         title    = title,
         subtitle = subtitle) +
    theme_fig() +
    theme(axis.title.y = if (show_y) element_text() else element_blank())
}

# ── Panel builder: with salience split (Studies 4, 5) ───────────────────────
make_acc_panel_sal <- function(df, title, subtitle = NULL,
                                show_y = TRUE, show_legend = FALSE,
                                x_breaks = c(-0.5, 0, 0.5),
                                x_labels = c("-0.5", "0", "0.5"),
                                y_floor = 0.33) {
  df <- dplyr::filter(df,
                      Direction %in% c("Faster", "Slower"),
                      !is.na(Accuracy), !is.na(Change), !is.na(Salience)) |>
    dplyr::mutate(Salience = factor(Salience, levels = c("High", "Low")))

  m <- tryCatch(
    lme4::glmer(Accuracy ~ poly(Change, 2) * Salience + (1 | id),
                data = df, family = binomial,
                control = glmerControl(optimizer = "bobyqa")),
    error = function(e)
      lme4::glmer(Accuracy ~ Change * Salience + (1 | id),
                  data = df, family = binomial,
                  control = glmerControl(optimizer = "bobyqa"))
  )

  pred <- ggeffects::ggpredict(m, terms = c("Change [all]", "Salience"))

  ggplot(pred, aes(x = x, y = predicted,
                   colour = group, fill = group)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                alpha = 0.20, colour = NA) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = sal_cols, labels = sal_labels,
                        name = "Salience") +
    scale_fill_manual(values = sal_cols, labels = sal_labels,
                      name = "Salience") +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = c(y_floor, 1),
                       breaks = c(0.33, 0.67, 1),
                       labels = c(".33", ".67", "1")) +
    labs(x        = "Breathing rate change",
         y        = if (show_y) "P(correct detection)" else NULL,
         title    = title,
         subtitle = subtitle) +
    theme_fig() +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      axis.title.y    = if (show_y) element_text() else element_blank()
    )
}

# ── Load data ────────────────────────────────────────────────────────────────
s1l <- readr::read_csv(file.path(DATA_DIR, "study1_long.csv")) |>
  dplyr::filter(Group == "TaskA")
s2l <- readr::read_csv(file.path(DATA_DIR, "study2_long.csv"))
s3l <- readr::read_csv(file.path(DATA_DIR, "study3_long.csv"))
s4l <- readr::read_csv(file.path(DATA_DIR, "study4_long.csv"))
s5l <- readr::read_csv(file.path(DATA_DIR, "study5_long.csv"))

# ── Build panels ─────────────────────────────────────────────────────────────
p1 <- make_acc_panel(s1l,
                     title    = "Study 1",
                     subtitle = "Online  |  N = 181",
                     show_y   = TRUE)

p2 <- make_acc_panel(s2l,
                     title    = "Study 2",
                     subtitle = "Online  |  N = 166",
                     show_y   = FALSE)

p3 <- make_acc_panel_sal(s3l,
                         title    = "Study 3 *",
                         subtitle = "Lab  |  N = 103  |  Fixed magnitudes",
                         show_y   = FALSE,
                         x_breaks = c(-0.65, 0, 0.65),
                         x_labels = c("-0.65", "0", "0.65"),
                         y_floor  = 0.10)

p4 <- make_acc_panel_sal(s4l,
                          title    = "Study 4",
                          subtitle = "Online  |  N = 131",
                          show_y   = TRUE)

p5 <- make_acc_panel_sal(s5l |>
                            dplyr::filter(Condition == "breath"),
                          title       = "Study 5",
                          subtitle    = "Lab  |  N = 206",
                          show_y      = FALSE,
                          show_legend = TRUE)

# ── Forest plot (Panel F) ────────────────────────────────────────────────────
# Shows partial r of Change² on Accuracy per study.
# Studies 4-5 split by Salience (High = orange, Low = blue) to show
# that high salience produces a steeper accuracy-magnitude curve.
# Studies 1, 2, 3: single estimate (no salience manipulation).
message("Computing Change² -> Accuracy partial r for forest plot...")

.changeq_row <- function(df, study_label, salience_label = "Overall", n_pp) {
  df <- df |>
    dplyr::filter(Direction %in% c("Faster", "Slower"),
                  !is.na(Accuracy), !is.na(Change)) |>
    dplyr::mutate(change_sq = Change^2)
  r <- .partial_r_glmer(df, "change_sq")
  dplyr::mutate(r, study = study_label, salience = salience_label, n = n_pp)
}

forest_df <- dplyr::bind_rows(
  # Studies without salience split
  .changeq_row(s1l, "Study 1",   n_pp = 181),
  .changeq_row(s2l, "Study 2",   n_pp = 166),
  # Study 3: single pooled estimate (fixed magnitudes inflate partial r;
  # salience split not comparable to staircase studies)
  .changeq_row(dplyr::filter(s3l, Direction %in% c("Faster","Slower")),
               "Study 3 *", n_pp = 103),
  # Studies 4-5 split by salience
  .changeq_row(dplyr::filter(s4l, Salience == "High"), "Study 4", "High", 131),
  .changeq_row(dplyr::filter(s4l, Salience == "Low"),  "Study 4", "Low",  131),
  .changeq_row(dplyr::filter(s5l, Condition == "breath", Salience == "High"),
               "Study 5", "High", 206),
  .changeq_row(dplyr::filter(s5l, Condition == "breath", Salience == "Low"),
               "Study 5", "Low",  206)
) |>
  dplyr::mutate(
    is_fixed  = study == "Study 3 *",
    salience  = factor(salience, levels = c("High", "Low", "Overall"))
  )

# RE pooled estimate (Fisher z, weighted by 1/se²)
pool_rows <- forest_df |>
  dplyr::filter(!is.na(r)) |>
  dplyr::mutate(se_r = sqrt((1 - r^2)^2 / (n - 2)),
                z    = atanh(r),
                wt   = 1 / se_r^2)
z_pool  <- sum(pool_rows$z * pool_rows$wt) / sum(pool_rows$wt)
se_pool <- sqrt(1 / sum(pool_rows$wt))
pooled_row <- tibble::tibble(
  study    = "Pooled", salience = factor("Overall", levels = c("High","Low","Overall")),
  r        = round(tanh(z_pool), 3),
  lower    = round(tanh(z_pool - 1.96 * se_pool), 3),
  upper    = round(tanh(z_pool + 1.96 * se_pool), 3),
  is_fixed = FALSE, n = NA_integer_
)
forest_df <- dplyr::bind_rows(forest_df, pooled_row)

display_order <- c("Study 1","Study 2","Study 3 *","Study 4","Study 5","Pooled")
forest_df <- forest_df |>
  dplyr::mutate(study = factor(study, levels = display_order))

sal_forest_cols <- c("High" = col_high, "Low" = col_low, "Overall" = "black")

p_forest <- ggplot(forest_df,
                   aes(x = r, y = study,
                       xmin = lower, xmax = upper,
                       colour = salience, shape = salience)) +
  scale_y_discrete(limits = rev(display_order),
                   expand  = expansion(add = 0.7)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = 1.5, colour = "grey40", linewidth = 0.5) +
  # Study-level error bars and points
  geom_errorbarh(
    data     = function(d) dplyr::filter(d, study != "Pooled"),
    height   = 0.18, linewidth = 0.6,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    data     = function(d) dplyr::filter(d, study != "Pooled", !is_fixed),
    size     = 3,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    data     = function(d) dplyr::filter(d, is_fixed),
    shape    = 17, size = 3,
    position = position_dodge(width = 0.6)
  ) +
  # Pooled
  geom_errorbarh(
    data     = function(d) dplyr::filter(d, study == "Pooled"),
    height   = 0.18, linewidth = 1.2,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    data     = function(d) dplyr::filter(d, study == "Pooled"),
    shape    = 18, size = 5,
    position = position_dodge(width = 0.6)
  ) +
  scale_colour_manual(values = sal_forest_cols,
                      labels = c("High" = "High salience",
                                 "Low"  = "Low salience",
                                 "Overall" = "No salience manipulation"),
                      name = NULL) +
  scale_shape_manual(values  = c("High" = 15, "Low" = 16, "Overall" = 15),
                     name    = NULL) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4),
                     labels = c("0", ".1", ".2", ".3", ".4")) +
  labs(
    x        = expression(paste("Partial ", italic(r))),
    y        = NULL,
    title    = "Change\u00b2 Predicts Accuracy",
    subtitle = "95% CI  |  RE meta-analysis"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 9, colour = "grey40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 10),
    legend.position    = "bottom",
    legend.text        = element_text(size = 9)
  )

# ── Combine ──────────────────────────────────────────────────────────────────
p_combined <- (p1 | p2 | p3) / (p4 | p5 | p_forest) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    caption    = paste0(
      "* Study 3 used fixed change magnitudes (+/-.20 to +/-.65); ",
      "partial r is inflated relative to staircase studies.\n",
      "Shaded regions: 95% CI.  Studies 3-5: High (orange) and Low (blue) salience.\n",
      "Panel F: partial r of Change\u00b2 on accuracy; Studies 4-5 split by salience; ",
      "Study 3 pooled (salience split not comparable due to fixed magnitudes)."
    ),
    theme = theme(
      plot.tag        = element_text(face = "bold", size = 12),
      plot.caption    = element_text(size = 8.5, colour = "grey40",
                                    hjust = 0, margin = margin(t = 6)),
      legend.position = "bottom"
    )
  )



print(p_combined)

ggsave(file.path(FIG_DIR, "fig_accuracy.pdf"),
       plot = p_combined, width = 10, height = 7, device = "pdf")
ggsave(file.path(FIG_DIR, "fig_accuracy.png"),
       plot = p_combined, width = 10, height = 7, dpi = 300)
cat("Saved: fig_accuracy\n")
