# ============================================================
# theme_bcat.R
# Persistent colour palette and ggplot theme for BCAT paper
#
# Source this file at the top of any script that produces figures.
# Conventions:
#   ORANGE (#ea801c) = "positive" condition:
#     Hit / Accurate / Detected / High salience / Faster / Acceleration
#   BLUE   (#1a80bb) = "negative" condition:
#     Miss / Inaccurate / Undetected / Low salience / Slower / Deceleration
# ============================================================

# ── Core palette ──────────────────────────────────────────────
BCAT_ORANGE <- "#ea801c"
BCAT_BLUE   <- "#1a80bb"

bcat_colours <- c(
  # Accuracy / Detection
  "Detected"     = BCAT_ORANGE,
  "Missed"       = BCAT_BLUE,
  "Hit"          = BCAT_ORANGE,
  "Miss"         = BCAT_BLUE,
  "1"            = BCAT_ORANGE,   # Accuracy == 1
  "0"            = BCAT_BLUE,     # Accuracy == 0
  # Salience
  "High"         = BCAT_ORANGE,
  "Low"          = BCAT_BLUE,
  "high"         = BCAT_ORANGE,
  "low"          = BCAT_BLUE,
  # Direction
  "Faster"       = BCAT_ORANGE,
  "Slower"       = BCAT_BLUE,
  "Acceleration" = BCAT_ORANGE,
  "Deceleration" = BCAT_BLUE
)

# ── Convenience scale functions ───────────────────────────────
# Drop these into any ggplot call like a normal scale_*

scale_colour_bcat <- function(name = NULL, ...) {
  ggplot2::scale_colour_manual(values = bcat_colours, name = name, ...)
}

scale_fill_bcat <- function(name = NULL, ...) {
  ggplot2::scale_fill_manual(values = bcat_colours, name = name, ...)
}

# ── Base theme ────────────────────────────────────────────────
theme_bcat <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      strip.text        = ggplot2::element_text(face = "bold"),
      legend.position   = "right",
      plot.title        = ggplot2::element_text(face = "bold", size = base_size),
      axis.title        = ggplot2::element_text(size = base_size - 1)
    )
}


# ── Study-level colour palette (for forest plots, multi-study figures) ────
BCAT_STUDY_COLOURS <- c(
  "Study1A" = "#4e9af1",
  "Study2"  = "#22c55e",
  "Study3"  = "#f59e0b",
  "Study4"  = "#ef4444",
  "Study5"  = "#8b5cf6"
)

# ── Figure save helper ────────────────────────────────────────
# Standard dimensions and PDF output for all BCAT paper figures.
# `path`: full file path including filename (no extension needed)
# `plot`: ggplot object; defaults to last plot
# `width`, `height`: inches
save_bcat_fig <- function(path, plot = ggplot2::last_plot(),
                           width = 7, height = 4.5) {
  path <- sub("\\.pdf$", "", path)
  ggplot2::ggsave(
    filename = paste0(path, ".pdf"),
    plot     = plot,
    width    = width,
    height   = height,
    device   = "pdf"
  )
  message("Saved: ", basename(paste0(path, ".pdf")))
  invisible(plot)
}
