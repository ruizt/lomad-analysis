# figure-theme.R -- one set of text sizes for every figure in the paper
#
# Sourced by mb-analysis/utils.R and simulations/simulation-results.R.
#
# Sizes are POINTS, chosen against the manuscript: body 11pt, captions 10pt
# (caption package, font=small). Figure text sits just under the caption.
#
# These land on the page as written only if the figure is rendered at the width
# it is placed at. \textwidth is 6.5in; a figure rendered at 10in and placed at
# 6.5in has all of its text scaled by 0.65.

PT <- list(
  title  = 9,     # axis titles
  text   = 7.5,   # tick labels
  strip  = 8,     # facet strips
  legend = 8,     # legend keys
  ltitle = 9,     # legend titles, panel titles
  annot  = 8      # in-panel text
)

# geom_text(), annotate("text") and geom_sf_text() take MILLIMETRES, not points.
# Everything drawn inside a panel must go through this.
ANNOT <- PT$annot / ggplot2::.pt      # 8pt -> 2.81mm

fig_sizes <- function() {
  ggplot2::theme(
    axis.title   = ggplot2::element_text(size = PT$title),
    axis.text    = ggplot2::element_text(size = PT$text),
    strip.text   = ggplot2::element_text(size = PT$strip),
    legend.text  = ggplot2::element_text(size = PT$legend),
    legend.title = ggplot2::element_text(size = PT$ltitle),
    plot.title   = ggplot2::element_text(size = PT$ltitle),
    plot.tag     = ggplot2::element_text(size = PT$ltitle, face = "bold"))
}

fig_theme <- function() {
  ggplot2::theme_minimal(base_size = PT$title) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.1, colour = "darkgray")) +
    fig_sizes()
}
