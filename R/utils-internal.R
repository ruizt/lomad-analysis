# Internal bootstrap progress display.
# Prints a single overwriting line: [===...] XX%  (b/B)
# Moves to a new line when b == B.
.boot_progress <- function(b, B, width = 40L) {
  filled <- round(b / B * width)
  bar    <- paste0(strrep("=", filled), strrep(" ", width - filled))
  cat(sprintf("\r  [%s] %3d%%  (%d/%d)", bar, round(b / B * 100), b, B))
  if (b == B) cat("\n")
  flush(stdout())
}
