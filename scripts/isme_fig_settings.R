# ISME Journal Figure Settings (Nature-style specifications)
# Save this script as: isme_fig_settings.R
# ------------------------------------------------------------
# Load this in your R script with: source("isme_fig_settings.R")
# ------------------------------------------------------------

# Recommended figure dimensions for single-column figures in ISME:
FIG_WIDTH  <- 85/25.4   # 85 mm = 3.35 in (single column)
FIG_HEIGHT <- 70/25.4   # 70 mm = ~2.75 in

# Optional for double-column figures:
FIG_WIDTH_DOUBLE <- 178/25.4  # 178 mm = 7.0 in (double column)

# Base font settings following Nature/ISME readability guidelines:
BASE_FONT <- "Arial"   # ISME prefers sans-serif fonts
BASE_SIZE <- 8         # ~8 pt text at final printed size

# Load required package
library(ggplot2)

# Set global ggplot2 theme for all figures
theme_set(
  theme_classic(base_size = BASE_SIZE, base_family = BASE_FONT) +
    theme(
      axis.text       = element_text(size = BASE_SIZE),
      axis.title      = element_text(size = BASE_SIZE + 1),
      legend.text     = element_text(size = BASE_SIZE),
      legend.title    = element_text(size = BASE_SIZE + 1),
      plot.title      = element_text(size = BASE_SIZE + 2, face = "bold"),
      strip.text      = element_text(size = BASE_SIZE + 1)
    )
)

# Helper function for saving ISME-ready figures
save_figure <- function(filename, plot, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = 300) {
  ggsave(filename, plot = plot, width = width, height = height, dpi = dpi, units = "in")
}
