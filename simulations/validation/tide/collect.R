## collect.R — assemble per-job .rds files and generate composite figure
##
## Run after fetch.sh has copied per-job .rds files locally.
##
## Usage (from the repo root):
##   Rscript simulations/validation/tide/collect.R
##
## Override the source directory:
##   RAW_DIR=/some/other/path Rscript simulations/validation/tide/collect.R

library(ggplot2)
library(patchwork)

RAW_DIR <- Sys.getenv("RAW_DIR", "simulations/validation/results/raw")
OUT_DIR <- "simulations/validation/results"

# ---- Load all results --------------------------------------------------------

files <- list.files(RAW_DIR, pattern = "\\.rds$", full.names = TRUE)

if (length(files) == 0) {
  stop("No .rds files found in: ", RAW_DIR,
       "\nCheck that the PVC is mounted and all jobs have completed.")
}

results <- lapply(files, readRDS)
names(results) <- vapply(results, function(r) r$experiment, character(1))

cat(sprintf("Loaded %d experiment files: %s\n",
            length(results), paste(names(results), collapse = ", ")))

# ---- Shared theme ------------------------------------------------------------

base_theme <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92")
  )

col_th  <- "firebrick"
col_emp <- "grey60"
col_oracle   <- "#0072B2"
col_pipeline <- "#D55E00"

# ==============================================================================
# Row (a): CLT QQ plots (oracle, faceted by s)
# ==============================================================================

clt_ids <- grep("^clt-", names(results), value = TRUE)

qq_list <- lapply(clt_ids, function(id) {
  r     <- results[[id]]
  R_vec <- r$R_vec[!is.na(r$R_vec)]
  Z     <- sqrt(r$s) * (R_vec - r$rho_oracle) / sqrt(r$V_oracle)
  nn    <- length(Z)

  data.frame(
    theoretical = qnorm(ppoints(nn)),
    empirical   = sort(Z),
    s_label     = sprintf("s = %d", r$s),
    s_num       = r$s
  )
})
qq_df <- do.call(rbind, qq_list)
qq_df$s_label <- reorder(qq_df$s_label, qq_df$s_num)

row_a <- ggplot(qq_df, aes(theoretical, empirical)) +
  geom_abline(slope = 1, intercept = 0, colour = col_th, linewidth = 0.5) +
  geom_point(colour = "black", alpha = 0.15, size = 0.9) +
  facet_wrap(~ s_label, nrow = 1) +
  labs(x = "Theoretical N(0,1)", y = expression("Empirical" ~ Z[t])) +
  base_theme

# ==============================================================================
# Row (b): Proposition 1 moment accuracy (oracle, s = 150)
# ==============================================================================

# ---- Left panel: rho --------------------------------------------------------

r_rho <- results[["rho-s150"]]
valid <- which(!is.na(r_rho$R_mean) & !is.na(r_rho$rho_th))

rho_df <- data.frame(
  t   = rep(valid, 2),
  rho = c(r_rho$rho_th[valid], r_rho$R_mean[valid]),
  type = rep(c("theoretical", "empirical"), each = length(valid))
)

p_rho <- ggplot() +
  geom_line(data = rho_df[rho_df$type == "empirical", ],
            aes(t, rho), colour = "grey30", linewidth = 0.5) +
  geom_line(data = rho_df[rho_df$type == "theoretical", ],
            aes(t, rho), colour = col_th, linewidth = 0.6) +
  # Dummy layer for legend
  geom_line(data = data.frame(
              x = c(NA, NA), y = c(NA, NA),
              label = factor(c("Empirical", "Theoretical"),
                             levels = c("Empirical", "Theoretical"))),
            aes(x, y, colour = label)) +
  scale_colour_manual(values = c("Empirical" = "grey30",
                                  "Theoretical" = col_th)) +
  labs(x = "Time", y = expression(rho[t])) +
  base_theme +
  theme(
    legend.position = c(1,1),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )

# ---- Right panel: V ----------------------------------------------------------

r_var   <- results[["var-s150"]]
V_emp   <- r_var$s * apply(r_var$R_mat, 2, var, na.rm = TRUE)
V_theory <- r_var$V_theory
valid_v  <- which(!is.na(V_emp) & !is.na(V_theory) & V_theory > 0)

v_df <- data.frame(
  V_theory = V_theory[valid_v],
  V_emp    = V_emp[valid_v]
)
rng <- range(c(v_df$V_theory, v_df$V_emp))

p_v <- ggplot(v_df, aes(V_theory, V_emp)) +
  geom_abline(slope = 1, intercept = 0, colour = col_th, linewidth = 0.6) +
  geom_point(colour = "black", size = 1.2, alpha = 0.4) +
  coord_equal(xlim = rng, ylim = rng) +
  labs(x = expression("Theoretical" ~ V[t]),
       y = expression(s %.% Var(R[t]))) +
  base_theme

row_b <- p_rho + p_v

# ==============================================================================
# Row (c): End-to-end pipeline validation (s = 150)
# ==============================================================================

r_e2e    <- results[["e2e-s150"]]
s_val    <- r_e2e$s
eval_pts <- r_e2e$eval_pts
n_pts    <- length(eval_pts)

# ---- Left panel: QQ at t = 1000 ---------------------------------------------

j_mid    <- which(eval_pts == 1000)
R_j      <- r_e2e$R_mat[, j_mid]
rho_or   <- r_e2e$rho_oracle[j_mid]
V_or     <- r_e2e$V_oracle[j_mid]
Z_oracle <- sqrt(s_val) * (R_j - rho_or) / sqrt(V_or)
Z_oracle <- Z_oracle[!is.na(Z_oracle)]
Z_pipe   <- r_e2e$Z_est_mat[, j_mid]
Z_pipe   <- Z_pipe[!is.na(Z_pipe)]

nn_e <- min(length(Z_oracle), length(Z_pipe))
qq_e <- data.frame(
  theoretical = rep(qnorm(ppoints(nn_e)), 2),
  empirical   = c(sort(Z_oracle[seq_len(nn_e)]), sort(Z_pipe[seq_len(nn_e)])),
  type        = rep(c("Oracle", "Pipeline"), each = nn_e)
)

p_qq <- ggplot(qq_e, aes(theoretical, empirical, colour = type)) +
  geom_abline(slope = 1, intercept = 0, colour = col_th, linetype = "dashed") +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_colour_manual(values = c(Oracle = col_oracle, Pipeline = col_pipeline)) +
  labs(x = "Theoretical N(0,1)", y = expression("Empirical" ~ Z[t])) +
  base_theme +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

# ---- Right panel: coverage ---------------------------------------------------

cov_rows <- vector("list", 2 * n_pts)
for (j in seq_len(n_pts)) {
  rho_j <- r_e2e$rho_oracle[j]
  V_j   <- r_e2e$V_oracle[j]

  Zo  <- sqrt(s_val) * (r_e2e$R_mat[, j] - rho_j) / sqrt(V_j)
  Zo  <- Zo[!is.na(Zo)]
  c_o <- mean(abs(Zo) < qnorm(0.975))
  n_o <- length(Zo)

  Zp  <- r_e2e$Z_est_mat[, j]
  Zp  <- Zp[!is.na(Zp)]
  c_p <- mean(abs(Zp) < qnorm(0.975))
  n_p <- length(Zp)

  cov_rows[[j]]         <- data.frame(
    t = eval_pts[j], cov = c_o,
    lo = c_o - 1.96 * sqrt(c_o * (1 - c_o) / n_o),
    hi = c_o + 1.96 * sqrt(c_o * (1 - c_o) / n_o),
    type = "Oracle"
  )
  cov_rows[[j + n_pts]] <- data.frame(
    t = eval_pts[j], cov = c_p,
    lo = c_p - 1.96 * sqrt(c_p * (1 - c_p) / n_p),
    hi = c_p + 1.96 * sqrt(c_p * (1 - c_p) / n_p),
    type = "Pipeline"
  )
}
cov_df <- do.call(rbind, cov_rows)

dodge <- position_dodge(width = 50)

p_cov <- ggplot(cov_df, aes(x = t, y = cov, colour = type)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = col_th) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 55,
                position = dodge, linewidth = 0.5) +
  geom_point(size = 1.5, position = dodge) +
  scale_colour_manual(values = c(Oracle = col_oracle, Pipeline = col_pipeline)) +
  labs(x = "Time", y = "95% coverage") +
  base_theme +
  theme(legend.position = "none")

row_c <- p_qq + p_cov

# ==============================================================================
# Composite figure
# ==============================================================================

# Use design layout so patchwork can align axes across rows.
# Row tags are added via labs(tag) on the first panel of each row.
row_a <- row_a + labs(tag = "a")
p_rho <- p_rho + labs(tag = "b")
p_qq  <- p_qq  + labs(tag = "c")

design <- "
AAAAAA
AAAAAA
BBBBCC
BBBBCC
DDDEEE
DDDEEE
DDDEEE
"

full_fig <- row_a + p_rho + p_v + p_qq + p_cov +
  plot_layout(design = design) +
  plot_annotation(
    theme = theme(plot.tag = element_text(size = 12, face = "bold"))
  )

ggsave(file.path(OUT_DIR, "validation-composite.png"),
       full_fig, width = 6, height = 6, dpi = 300)
# ggsave(file.path(OUT_DIR, "validation-composite.pdf"),
#        full_fig, width = 10, height = 9)

cat("Saved validation-composite.png and .pdf\n")
cat("\nDone.\n")
