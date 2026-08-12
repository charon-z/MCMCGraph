# fig_style.R --------------------------------------------------------------
# Shared figure design system for BPFC manuscript figures, so simulation
# and real-data (mouse, LINCS) figures look identical and journal-grade.
#
# Source this AFTER library(ggplot2). Provides:
#   MCG_FONT, MCG_BASE          font family + base size
#   MCG_PAL_IND                 the two paired states (Individual 1 / 2): orange / blue
#   MCG_PAL_QUAL                Okabe-Ito colourblind-safe qualitative palette (clusters/methods)
#   mcg_theme()                 the shared ggplot theme
#   mcg_open_pdf() / mcg_open_png()   macOS quartz devices for real-Arial, 300-dpi output
#   mcg_save(name, plot/draw, w, h)   save both .pdf and .png with the shared devices

MCG_FONT <- "Arial"
MCG_BASE <- 26

# Okabe-Ito (Nature/Science-standard colourblind-safe qualitative palette)
MCG_OKABE <- c("#0072B2","#E69F00","#009E73","#CC79A7","#D55E00","#56B4E9","#F0E442","#000000")
MCG_PAL_QUAL <- MCG_OKABE
# Two paired states use a fixed convention everywhere: state 1 = orange, state 2 = blue
# (vermillion orange reads cleanly over dense faint trajectories; both are Okabe-Ito)
MCG_PAL_IND  <- c("#D55E00", "#0072B2")

# Canonical 15-module palette (Dark2 + Paired). SINGLE SOURCE OF TRUTH: every
# figure that colours by module (multilayer network, GO enrichment, ...) must
# take colours from here, keyed M1..M15, so module identity == colour everywhere.
MCG_MOD_PAL <- c(M1="#1B9E77", M2="#D95F02", M3="#7570B3", M4="#E7298A", M5="#66A61E",
                 M6="#E6AB02", M7="#A6761D", M8="#666666", M9="#1F78B4", M10="#B2DF8A",
                 M11="#FB9A99", M12="#FDBF6F", M13="#CAB2D6", M14="#FF7F00", M15="#6A3D9A")

mcg_theme <- function(base_size = MCG_BASE, legend = "none") {
  ggplot2::theme_bw(base_size = base_size, base_family = MCG_FONT) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.3, colour = "grey88"),
      panel.border       = ggplot2::element_rect(linewidth = 0.7, colour = "grey25"),
      axis.ticks         = ggplot2::element_line(linewidth = 0.5, colour = "grey25"),
      axis.title         = ggplot2::element_text(),
      plot.title         = ggplot2::element_text(hjust = 0),
      plot.subtitle      = ggplot2::element_text(colour = "grey30"),
      strip.background   = ggplot2::element_rect(fill = "grey92", colour = NA),
      strip.text         = ggplot2::element_text(),
      legend.position    = legend,
      legend.key         = ggplot2::element_blank(),
      plot.margin        = ggplot2::margin(8, 10, 6, 6)
    )
}

# macOS quartz devices => real Arial (CoreText) + sharp 300 dpi. Do NOT use cairo_pdf
# (it silently substitutes Arial -> Hiragino on this fontconfig).
mcg_save <- function(name, obj, w = 9, h = 5, dir = ".") {
  is_fun <- is.function(obj)
  grDevices::quartz(file = file.path(dir, paste0(name, ".pdf")), type = "pdf",
                    width = w, height = h, family = MCG_FONT)
  if (is_fun) obj() else print(obj); grDevices::dev.off()
  grDevices::png(file.path(dir, paste0(name, ".png")), width = w, height = h,
                 units = "in", res = 400, type = "quartz", family = MCG_FONT)
  if (is_fun) obj() else print(obj); grDevices::dev.off()
}

# base-graphics defaults matching the theme (call inside a draw() before plotting)
mcg_par <- function(...) {
  graphics::par(family = MCG_FONT, cex.main = 1.7, cex.lab = 1.6, cex.axis = 1.35,
                font.main = 2, font.lab = 2, col.axis = "grey25", fg = "grey25", ...)
}
