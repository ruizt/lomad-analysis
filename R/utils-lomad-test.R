# Internal test implementation for lomad_test()
#
# .lomad_test_clt() — pointwise CLT test (paper method)


# Pointwise CLT test: the paper method.
.lomad_test_clt <- function(fit, alpha) {
  if (alpha <= 0 || alpha >= 1)
    stop("`alpha` must be in (0, 1).")

  n         <- fit$inputs$n
  s         <- fit$inputs$s
  valid_idx <- fit$valid_idx
  m         <- length(valid_idx)

  if (m == 0L) stop("No valid time points in fit.")

  R   <- fit$R
  rho <- fit$rho
  V   <- fit$V

  Z      <- rep(NA_real_, n)
  p_raw  <- rep(NA_real_, n)

  se                  <- sqrt(V[valid_idx] / s)
  Z[valid_idx]        <- (R[valid_idx] - rho[valid_idx]) / se
  p_raw[valid_idx]    <- stats::pnorm(Z[valid_idx])

  alpha_eff           <- alpha / sum(1 / seq_len(m))
  p_adj               <- rep(NA_real_, n)
  p_adj[valid_idx]    <- pmin(p_raw[valid_idx] * sum(1 / seq_len(m)), 1)

  rejected            <- rep(NA, n)
  rejected[valid_idx] <- p_raw[valid_idx] <= alpha_eff

  n_rejected <- sum(rejected, na.rm = TRUE)
  message(sprintf("Rejected %d / %d time points at BY-FDR = %.2f (alpha_eff = %.4f)",
                  n_rejected, m, alpha, alpha_eff))

  list(
    Z         = Z,
    p_values  = p_raw,
    p_adj     = p_adj,
    rejected  = rejected,
    alpha_eff = alpha_eff,
    inputs    = list(alpha = alpha, s = s)
  )
}
