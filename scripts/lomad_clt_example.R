library(lomad)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# ---- Parameters -------------------------------------------------------------

n      <- 1000
d      <- 1.5     # L2 separation between trend coefficient vectors
phi    <- 0.5     # AR(1) autocorrelation
snr    <- 1.5     # target signal-to-noise ratio (passed to sim_noise_pair)
h      <- NULL    # MA smoothing window  (NULL = auto: max(5, floor(n/200)))
s      <- NULL    # rolling correlation  (NULL = auto: see lomad_fit)
seed   <- 3817

# ---- Generate unstructured trends -------------------------------------------
#
# "Unstructured" = two Fourier-basis series separated by a fixed L2 distance d,
# with no time-varying coupling weight w_t.  The series are globally distinct
# but there is no imposed temporal structure to where they decouple.

trends <- sim_trends(n = n, d = d, method = "dist", seed = seed)

# ---- Add AR(1) noise --------------------------------------------------------

# sim_noise_pair needs a concrete h for SNR calibration; resolve the auto value
# here so it matches what lomad_fit will use internally.
h_eff <- if (is.null(h)) max(5L, floor(n / 200L)) else as.integer(h)

sim <- sim_noise_pair(
  trends        = trends,
  h             = h_eff,
  lambda_target = snr,
  ar.coefs      = phi,
  seed          = seed
)

cat(sprintf(
  "SNR  series 1: %.2f   series 2: %.2f\n",
  sim$noise$series1$mean_snr,
  sim$noise$series2$mean_snr
))

# ---- Fit + test -------------------------------------------------------------

fit <- lomad_fit(sim$y1, sim$y2, h = h, s = s, max_pq = 3L)
tst <- lomad_test(fit, alpha = 0.05)

cat(sprintf(
  "ARMA fit: (%d,%d) / (%d,%d)   innovation sigma2: %.3f / %.3f\n",
  fit$noise$series1$order[1], fit$noise$series1$order[3],
  fit$noise$series2$order[1], fit$noise$series2$order[3],
  fit$noise$series1$sigma2,   fit$noise$series2$sigma2
))

vi <- fit$valid_idx
cat(sprintf("Rejected %d / %d valid time points (alpha = 0.05, BY-FDR)\n",
            sum(tst$rejected[vi]), length(vi)))

# ---- Plot -------------------------------------------------------------------

t_idx <- seq_len(n)

df_series <- tibble(
  t    = t_idx,
  y1   = sim$y1,  y2   = sim$y2,
  sm1  = fit$ma1, sm2  = fit$ma2
)

p_obs <- ggplot(df_series, aes(x = t)) +
  geom_line(aes(y = y1), colour = "#BBBBBB", linewidth = 0.4) +
  geom_line(aes(y = y2), colour = "#BBBBBB", linewidth = 0.4) +
  geom_line(aes(y = sm1, colour = "Series 1"), linewidth = 0.8) +
  geom_line(aes(y = sm2, colour = "Series 2"), linewidth = 0.8) +
  scale_colour_manual(values = c("Series 1" = "#0072B2", "Series 2" = "#D55E00")) +
  labs(x = NULL, y = "Value", colour = NULL,
       title = sprintf("Unstructured trends  (d = %.1f, AR(1) \u03d5 = %.1f)", d, phi)) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", plot.title = element_text(size = 11))

df_corr <- tibble(
  t   = t_idx,
  R   = fit$R,
  rho = fit$rho
) |> filter(t %in% vi)

p_corr <- ggplot(df_corr, aes(x = t)) +
  geom_line(aes(y = R,   colour = "R_t"),      linewidth = 0.7) +
  geom_line(aes(y = rho, colour = "\u03c1_t"), linewidth = 0.8, linetype = "dashed") +
  scale_colour_manual(values = c("R_t" = "#009E73", "\u03c1_t" = "#CC79A7")) +
  labs(x = NULL, y = "Correlation", colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

df_z <- tibble(
  t        = t_idx,
  Z        = tst$Z,
  rejected = tst$rejected
) |> filter(t %in% vi)

z_thresh <- qnorm(tst$alpha_eff)

p_z <- ggplot(df_z, aes(x = t)) +
  geom_hline(yintercept = 0,        colour = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = z_thresh, colour = "#E69F00", linewidth = 0.6,
             linetype = "dashed") +
  geom_line(aes(y = Z), colour = "#56B4E9", linewidth = 0.6) +
  geom_point(data = filter(df_z, isTRUE(rejected)),
             aes(y = Z), colour = "#D55E00", size = 1.5) +
  annotate("text", x = max(vi) * 0.02, y = z_thresh + 0.15,
           label = "BY threshold", hjust = 0, size = 3, colour = "#E69F00") +
  labs(x = "Time", y = "Z-statistic") +
  theme_minimal(base_size = 11)

fig <- p_obs / p_corr / p_z +
  plot_layout(heights = c(2, 1.5, 1.5)) &
  theme(panel.grid.minor = element_blank())

print(fig)

ggsave("scripts/_img/fig_clt_example.png",
       fig, width = 7, height = 6, dpi = 150)
