#' Pointwise CLT test for local correlation decoupling
#'
#' Given the output of [lomad_fit_clt()], tests at each valid time point
#' whether the local sample correlation \eqn{R_t} is significantly below the
#' estimated local population correlation \eqn{\hat\rho_t}. The test statistic
#' is
#' \deqn{Z_t = \frac{R_t - \hat\rho_t}{\sqrt{\hat V_t / s}}}
#' which is asymptotically \eqn{N(0, 1)} under local similarity. One-sided
#' p-values \eqn{P(Z \leq Z_t)} are computed and adjusted for multiple
#' comparisons using the Benjamini--Yekutieli (BY) procedure, which controls
#' FDR under arbitrary dependence.
#'
#' @param fit List returned by [lomad_fit_clt()].
#' @param alpha Numeric in \eqn{(0, 1)}. FDR level. Default 0.05.
#'
#' @return A named list:
#'   \describe{
#'     \item{Z}{Numeric vector (length `n`). Pointwise z-statistics; `NA`
#'       outside `valid_idx`.}
#'     \item{p_values}{Numeric vector (length `n`). Raw one-sided p-values;
#'       `NA` outside `valid_idx`.}
#'     \item{p_adj}{Numeric vector (length `n`). BY-adjusted p-values; `NA`
#'       outside `valid_idx`.}
#'     \item{rejected}{Logical vector (length `n`). `TRUE` where the local
#'       similarity hypothesis is rejected at level `alpha` after BY
#'       adjustment.}
#'     \item{alpha_eff}{Numeric. Effective BY p-value threshold.}
#'     \item{inputs}{List of `alpha` and `s` from the fit.}
#'   }
#'
#' @seealso [lomad_fit_clt()]
#'
#' @export
lomad_test_clt <- function(fit, alpha = 0.05) {
  if (alpha <= 0 || alpha >= 1)
    stop("`alpha` must be in (0, 1).")

  n         <- fit$inputs$n
  s         <- fit$inputs$s
  valid_idx <- fit$valid_idx
  m         <- length(valid_idx)

  if (m == 0L) stop("No valid time points in fit. Check `lomad_fit_clt()`.")

  R   <- fit$R
  rho <- fit$rho
  V   <- fit$V

  # --- Pointwise z-statistics and p-values ---
  Z      <- rep(NA_real_, n)
  p_raw  <- rep(NA_real_, n)

  se                  <- sqrt(V[valid_idx] / s)
  Z[valid_idx]        <- (R[valid_idx] - rho[valid_idx]) / se
  p_raw[valid_idx]    <- stats::pnorm(Z[valid_idx])  # lower tail: H1: R_t < rho_t

  # --- BY adjustment (controls FDR under arbitrary dependence) ---
  # alpha_eff = alpha / sum(1/1 + 1/2 + ... + 1/m)
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
