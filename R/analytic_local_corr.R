#' Analytic CLT inference for the fraction of time in the low-correlation state
#'
#' Provides a closed-form CLT-based p-value for the fraction of time in the
#' low-correlation state using a 2-state Markov chain approximation derived
#' from [lomad_fit()]. No simulation is required.
#'
#' @param fit_obj List returned by [lomad_fit()].
#' @param T_eff Integer or NULL. Effective chain length to use in the CLT. If
#'   `NULL`, uses `length(fit_obj$valid_idx)`.
#'
#' @return A list containing:
#'   \describe{
#'     \item{frac_state}{List with asymptotic mean (`pi1`), observed value
#'       (`obs`), CLT variance (`var_CLT`), standard error (`se_CLT`),
#'       z-statistic (`z`), one-sided p-value (`p_right`), and `T_eff`.}
#'     \item{transitions}{Markov transition parameters used: `pi1`, `pi0`,
#'       `pi11`, `pi01`, `pi00`.}
#'   }
#'
#' @note Currently provides a CLT-based p-value only for the fraction of time
#'   in the low-correlation state. Asymptotic expectations for entry rate and
#'   run length are given by `fit_obj$expected_asymptotic`.
#'
#' @export
analytic_local_corr <- function(fit_obj,
                                T_eff = NULL) {

  message("Starting analytic Markov CLT inference (frac_state only)...")

  observed  <- fit_obj$observed
  expectedA <- fit_obj$expected_asymptotic

  pi1_th  <- expectedA$frac_state
  L_th    <- expectedA$mean_run_length
  pi11_th <- if (!is.na(L_th) && L_th > 1) 1 - 1 / L_th else NA_real_
  pi0_th  <- if (!is.na(pi1_th)) 1 - pi1_th else NA_real_
  pi01_th <- if (!is.na(pi1_th) && !is.na(pi11_th) &&
                 pi1_th > 0 && pi1_th < 1) {
    pi1_th * (1 - pi11_th) / (1 - pi1_th)
  } else NA_real_
  pi00_th <- if (!is.na(pi01_th)) 1 - pi01_th else NA_real_

  if (is.null(T_eff)) {
    T_eff <- length(fit_obj$valid_idx)
  }

  ## CLT for sample mean of a 2-state Markov chain:
  ##   Var(I_t) = pi1 (1 - pi1)
  ##   theta = pi11 + pi00 - 1  (lag-1 autocorr of I_t)
  ##   Var(pi1_hat) ~ Var(I_t) * (1 + theta) / (1 - theta) / T_eff
  var_CLT <- NA_real_
  se_CLT  <- NA_real_
  z       <- NA_real_
  p_right <- NA_real_

  if (is.finite(pi1_th) && pi1_th > 0 && pi1_th < 1 &&
      is.finite(pi11_th) && is.finite(pi00_th) &&
      T_eff > 0) {

    var_I <- pi1_th * (1 - pi1_th)
    theta <- pi11_th + pi00_th - 1  # corr(I_t, I_{t+1})

    if (is.finite(theta) && abs(theta) < 1) {
      var_CLT <- var_I * (1 + theta) / (1 - theta) / T_eff
      if (var_CLT > 0) {
        se_CLT <- sqrt(var_CLT)
        z      <- (observed$frac_state - pi1_th) / se_CLT
        p_right <- 1 - stats::pnorm(z)
      }
    }
  }

  transitions <- list(
    pi1  = pi1_th,
    pi0  = pi0_th,
    pi11 = pi11_th,
    pi01 = pi01_th,
    pi00 = pi00_th
  )

  frac_state <- list(
    pi1     = pi1_th,
    obs     = observed$frac_state,
    var_CLT = var_CLT,
    se_CLT  = se_CLT,
    z       = z,
    p_right = p_right,
    T_eff   = T_eff
  )

  message("Analytic Markov CLT inference complete.")

  list(
    frac_state  = frac_state,
    transitions = transitions
  )
}
