# Combined replication figure
# Layout: 3 columns x 2 rows
#   Row 1: Study 1 | Study 2 | Study 3*
#   Row 2: Study 4 | Study 5 | Meta-analytic forest plot

# Set Up ---------

ctrl <- lmerControl(optimizer = "bobyqa")

# --- Colour theme -----------------------------------------------------------
# UPDATE these two hex values to match your theme
col_detected <- "#E69F00"   # orange — detected / hit trials
col_missed   <- "#0072B2"   # blue   — missed trials

acc_cols   <- c("0" = col_missed, "1" = col_detected)
acc_labels <- c("0" = "Missed",   "1" = "Detected")

# --- Helper: fit model and return ggplot panel ------------------------------

make_panel <- function(df, title, subtitle = NULL,
                       show_y = TRUE, show_legend = FALSE,
                       x_breaks = c(-0.5, 0, 0.5),
                       x_labels = c("-0.5", "0", "0.5")) {

  df$Arousal_z <- as.numeric(scale(df$Arousal))

  lm.full <- tryCatch(
    lmer(Arousal_z ~ Accuracy * poly(Change, 2) + (Change | id),
         data = df, control = ctrl),
    error = function(e)
      lmer(Arousal_z ~ Accuracy * poly(Change, 2) + (1 | id),
           data = df, control = ctrl)
  )
  if (isSingular(lm.full))
    lm.full <- lmer(Arousal_z ~ Accuracy * poly(Change, 2) + (1 | id),
                    data = df, control = ctrl)

  pred <- ggpredict(lm.full, terms = c("Change [all]", "Accuracy [0, 1]"))

  plot(pred) +
    scale_colour_manual(values = acc_cols, labels = acc_labels,
                        name = "Detection") +
    scale_fill_manual(values = acc_cols, labels = acc_labels,
                      name = "Detection") +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    labs(
      x        = "Breathing rate change",
      y        = if (show_y) "Predicted arousal (standardised)" else NULL,
      title    = title,
      subtitle = subtitle
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 11),
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      legend.position  = if (show_legend) "bottom" else "none",
      axis.title.y     = if (show_y) element_text() else element_blank(),
      panel.grid.minor = element_blank()
    )
}

# --- Load data --------------------------------------------------------------

s1l <- read.csv(file.path(DATA_DIR, "study1_long.csv"))
s2l <- read.csv(file.path(DATA_DIR, "study2_long.csv"))
s3l <- read.csv(file.path(DATA_DIR, "study3_long.csv"))
s4l <- read.csv(file.path(DATA_DIR, "study4_long.csv"))
s5l <- read.csv(file.path(DATA_DIR, "study5_long.csv"))

# --- Build study panels -----------------------------------------------------

p1 <- make_panel(s1l,
                 title    = "Study 1",
                 subtitle = "Online  |  N = 181",
                 show_y   = TRUE)

p2 <- make_panel(s2l,
                 title    = "Study 2",
                 subtitle = "Online  |  N = 166",
                 show_y   = FALSE)

# Study 3: fixed magnitudes — x-axis anchored at ±0.65
p3 <- make_panel(s3l,
                 title    = "Study 3 *",
                 subtitle = "Lab  |  N = 103  |  Fixed magnitudes",
                 show_y   = FALSE,
                 x_breaks = c(-0.65, 0, 0.65),
                 x_labels = c("-0.65", "0", "0.65"))

p4 <- make_panel(s4l,
                 title    = "Study 4",
                 subtitle = "Online  |  N = 131",   # corrected from Lab
                 show_y   = TRUE)

p5 <- make_panel(s5l,
                 title       = "Study 5",
                 subtitle    = "Lab  |  N = 206",
                 show_y      = FALSE,
                 show_legend = TRUE)

# --- Forest plot (Panel F) --------------------------------------------------
# display_order: studies first, Pooled last
# scale_y_discrete(limits = rev(display_order)) puts Pooled at bottom

forest_df <- data.frame(
  label     = c("Study 1", "Study 2", "Study 3 *",
                "Study 4", "Study 5", "Pooled"),
  r_val     = c(-0.110, -0.109, -0.184, -0.076, -0.058, -0.107),
  r_lower   = c(-0.183, -0.178, -0.202, -0.145, -0.093, -0.151),
  r_upper   = c(-0.037, -0.040, -0.166, -0.007, -0.023, -0.064),
  is_pooled = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  is_fixed  = c(FALSE, FALSE, TRUE,  FALSE, FALSE, FALSE)
)

display_order <- c("Study 1", "Study 2", "Study 3 *",
                   "Study 4", "Study 5", "Pooled")

p_forest <- ggplot(forest_df,
                   aes(x = r_val, y = label,
                       xmin = r_lower, xmax = r_upper)) +
  scale_y_discrete(
    limits = rev(display_order),
    expand = expansion(add = 0.7)   # vertical padding to match panel height
  ) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = 1.5, colour = "grey40", linewidth = 0.5) +
  # Individual study CIs — black
  geom_errorbarh(
    data   = function(d) dplyr::filter(d, !is_pooled),
    height = 0.2, linewidth = 0.7, colour = "black"
  ) +
  # Staircase study points — black squares
  geom_point(
    data  = function(d) dplyr::filter(d, !is_pooled, !is_fixed),
    shape = 15, size = 3.5, colour = "black"
  ) +
  # Fixed-magnitude study point — black triangle
  geom_point(
    data  = function(d) dplyr::filter(d, is_fixed),
    shape = 17, size = 3.5, colour = "black"
  ) +
  # Pooled CI — blue (missed colour)
  geom_errorbarh(
    data   = function(d) dplyr::filter(d, is_pooled),
    height = 0.2, linewidth = 1, colour = col_missed
  ) +
  # Pooled point — blue diamond
  geom_point(
    data  = function(d) dplyr::filter(d, is_pooled),
    shape = 18, size = 5, colour = col_missed
  ) +
  scale_x_continuous(
    limits = c(-0.24, 0.04),
    breaks = c(-0.2, -0.1, 0.0),
    labels = c("-0.2", "-0.1", "0.0")
  ) +
  labs(
    x        = expression(paste("Partial ", italic(r))),
    y        = NULL,
    title    = "H4B: Change x Accuracy",
    subtitle = "95% CI  |  RE meta-analysis"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 9, colour = "grey40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 10)
  )

# --- Combine ----------------------------------------------------------------

p_combined <- (p1 | p2 | p3) / (p4 | p5 | p_forest) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    caption    = paste0(
      "* Study 3 used fixed change magnitudes (+/-.20 to +/-.65); ",
      "partial r is inflated relative to staircase studies (see Supplementary S2.1).\n",
      "Arousal z-scored within study. Shaded regions: 95% CI."
    ),
    theme = theme(
      plot.tag        = element_text(face = "bold", size = 12),
      plot.caption    = element_text(size = 8.5, colour = "grey40",
                                     hjust = 0, margin = margin(t = 6)),
      legend.position = "bottom"
    )
  )

print(p_combined)

# --- Save -------------------------------------------------------------------



ggsave(file.path(FIG_DIR, "fig_arousal.pdf"),
       plot = p_combined, width = 10, height = 7, device = "pdf")

ggsave(file.path(FIG_DIR, "fig_arousal.png"),
       plot = p_combined, width = 10, height = 7, dpi = 300)

cat("Figure saved to", FIG_DIR, "\n")
