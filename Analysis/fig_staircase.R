# fig_staircase.R
# Staircase threshold level progression across studies
# Layout: 2 columns × 2 rows
#   Row 1: Study 1A | Study 2
#   Row 2: Study 4  | Study 5
#
# Studies 1A and 2: single 3AFC staircase — one line
# Studies 4 and 5:  four independent staircases —
#   2 Salience (High / Low) × 2 Direction (Faster / Slower)
#
# Y axis: mean staircase level (|Change|) per sequential trial step
#          averaged across participants, ± 1 SE ribbon
# X axis: sequential trial number within each staircase condition
#
# Colour:  Salience — High = orange, Low = blue (matches fig_accuracy.R)
# Linetype: Direction — Faster = solid, Slower = dashed

# Set Up ---------

# ── Colour / linetype palette ────────────────────────────────────────────────
col_high <- "#E69F00"   # High salience (orange)
col_low  <- "#0072B2"   # Low  salience (blue)

sal_cols    <- c("High" = col_high, "Low" = col_low)
sal_labels  <- c("High" = "High salience", "Low" = "Low salience")
dir_types   <- c("Faster" = "solid", "Slower" = "dashed")
dir_labels  <- c("Faster" = "Acceleration", "Slower" = "Deceleration")

# ── Shared theme ─────────────────────────────────────────────────────────────
theme_fig <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 11),
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      panel.grid.minor = element_blank(),
      legend.position  = "none"
    )
}

# ── Prep: compute within-condition trial number per participant ────────────────
# Returns per-participant trial-level data (not aggregated).
# geom_smooth handles the averaging and smoothing.
prep_staircase <- function(df, group_cols = character(0), max_trials = Inf) {
  df |>
    dplyr::filter(!is.na(abs_level)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c("id", group_cols)))) |>
    dplyr::arrange(trial_order, .by_group = TRUE) |>
    dplyr::mutate(trial_n = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::filter(trial_n <= max_trials)
}

# ── Load and prepare data ─────────────────────────────────────────────────────

# Study 1A: single staircase, no salience / direction split
# Use abs(Change) as staircase level; filter to change trials only
s1_raw <- readr::read_csv(file.path(DATA_DIR, "study1_long.csv")) |>
  dplyr::filter(Group == "TaskA",
                Direction %in% c("Faster", "Slower")) |>
  dplyr::mutate(abs_level   = abs(Change),
                trial_order = TrialNum)

s1_traj <- prep_staircase(s1_raw)

# Study 2: single staircase
s2_raw <- readr::read_csv(file.path(DATA_DIR, "study2_long.csv")) |>
  dplyr::filter(Direction %in% c("Faster", "Slower")) |>
  dplyr::mutate(abs_level   = abs(Change),
                trial_order = TrialNum)

s2_traj <- prep_staircase(s2_raw)

# Study 4: column is capitalised as "Level"
s4_df  <- readr::read_csv(file.path(DATA_DIR, "study4_long.csv"))
s4_raw <- s4_df |>
  dplyr::filter(Direction %in% c("Faster", "Slower"), !is.na(Salience)) |>
  dplyr::mutate(
    abs_level   = Level,
    trial_order = Trial,
    Salience    = factor(Salience, levels = c("High", "Low")),
    Direction   = factor(Direction, levels = c("Faster", "Slower"))
  ) |>
  dplyr::filter(!is.na(abs_level))

s4_traj <- prep_staircase(s4_raw, group_cols = c("Salience", "Direction"),
                           max_trials = 10)

# Study 5: breath condition only, restricted to each participant's FIRST
# breath session to avoid double-counting.
# Breath-first group: first breath session = ses1
# Visual-first group: first breath session = ses2
# Without this filter, breath-first participants contribute 20 trials
# (ses1 + ses2) per condition, pushing trial_n beyond 10.
s5_df  <- readr::read_csv(file.path(DATA_DIR, "study5_long.csv"))
s5_raw <- s5_df |>
  dplyr::filter(Condition == "breath",
                Direction %in% c("Faster", "Slower"),
                !is.na(Salience)) |>
  dplyr::group_by(id) |>
  dplyr::mutate(first_breath_ses = dplyr::first(sort(unique(ses)))) |>
  dplyr::ungroup() |>
  dplyr::filter(ses == first_breath_ses) |>
  dplyr::mutate(
    abs_level   = level,
    trial_order = Trial,
    Salience    = factor(Salience, levels = c("High", "Low")),
    Direction   = factor(Direction, levels = c("Faster", "Slower"))
  ) |>
  dplyr::filter(!is.na(abs_level))

s5_traj <- prep_staircase(s5_raw, group_cols = c("Salience", "Direction"),
                           max_trials = 10)

# ── Panel builders ────────────────────────────────────────────────────────────
# geom_smooth(method = "loess") smooths over the per-participant data.
# span = 1.0 for Studies 4-5 (only 10 trials; higher span = smoother).
# span = 0.75 for Studies 1A-2 (25+ trials; default loess span).

# Single-line panel (Studies 1A and 2) — blue (gradual / low-salience staircase)
make_staircase_single <- function(df, title, subtitle = NULL,
                                   show_y = TRUE, col = col_low,
                                   span = 0.75) {
  ggplot(df, aes(x = trial_n, y = abs_level)) +
    geom_smooth(method  = "loess", formula = y ~ x,
                span    = span,
                colour  = col, fill = col, alpha = 0.20,
                linewidth = 1, se = TRUE) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    labs(
      x     = "Trial",
      y     = if (show_y) "Staircase level (|\u0394rate|)" else NULL,
      title = title, subtitle = subtitle
    ) +
    theme_fig() +
    theme(axis.title.y = if (show_y) element_text() else element_blank())
}

# Four-line panel (Studies 4 and 5)
make_staircase_four <- function(df, title, subtitle = NULL,
                                 show_y = TRUE, show_legend = FALSE,
                                 span = 1.0) {
  ggplot(df,
         aes(x        = trial_n,
             y        = abs_level,
             colour   = Salience,
             fill     = Salience,
             linetype = Direction)) +
    geom_smooth(method    = "loess", formula = y ~ x,
                span      = span,
                alpha     = 0.15, linewidth = 1, se = TRUE) +
    scale_colour_manual(values = sal_cols, labels = sal_labels,
                        name = "Salience") +
    scale_fill_manual(values = sal_cols, labels = sal_labels,
                      name = "Salience") +
    scale_linetype_manual(values = dir_types, labels = dir_labels,
                          name = "Direction") +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    labs(
      x     = "Trial",
      y     = if (show_y) "Staircase level (|\u0394rate|)" else NULL,
      title = title, subtitle = subtitle
    ) +
    theme_fig() +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      axis.title.y    = if (show_y) element_text() else element_blank()
    )
}

# ── Build panels ──────────────────────────────────────────────────────────────
p1 <- make_staircase_single(
  s1_traj,
  title    = "Study 1A",
  subtitle = "Online  |  N = 181  |  Single staircase",
  show_y   = TRUE
)

p2 <- make_staircase_single(
  s2_traj,
  title    = "Study 2",
  subtitle = "Online  |  N = 166  |  Single staircase",
  show_y   = FALSE
)

p4 <- make_staircase_four(
  s4_traj,
  title    = "Study 4",
  subtitle = "Online  |  N = 131  |  High / Low salience \u00d7 Faster / Slower",
  show_y   = TRUE
)

p5 <- make_staircase_four(
  s5_traj,
  title       = "Study 5",
  subtitle    = "Lab  |  N = 206  |  Breath condition",
  show_y      = FALSE,
  show_legend = TRUE
)

# ── Combine ───────────────────────────────────────────────────────────────────



p_combined <- (p1 | p2) / (p4 | p5) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    caption    = paste0(
      "Studies 4-5: Quest staircase level (threshold estimate per trial), ",
      "averaged across participants (\u00b11 SE).\n",
      "Studies 1A-2: |Change| (abs. magnitude of signed change).\n",
      "Studies 4-5: High salience (orange) and Low salience (blue); ",
      "Acceleration (solid) and Deceleration (dashed)."
    ),
    theme = theme(
      plot.tag        = element_text(face = "bold", size = 12),
      plot.caption    = element_text(size = 8.5, colour = "grey40",
                                     hjust = 0, margin = margin(t = 6)),
      legend.position = "bottom"
    )
  )

print(p_combined)

ggsave(file.path(FIG_DIR, "fig_staircase.pdf"),
       plot = p_combined, width = 9, height = 7, device = "pdf")
ggsave(file.path(FIG_DIR, "fig_staircase.png"),
       plot = p_combined, width = 9, height = 7, dpi = 300)
cat("Saved: fig_staircase\n")
