## analytic_local_corr.R
##
## Analytic (CLT-based) inference for the fraction of time in the
## low-correlation state, using the 2-state Markov chain approximation.
##
## INPUTS:
##   fit_obj : list returned by fit_null_model_local_corr()
##   T_eff   : effective chain length to use in CLT; if NULL, uses
##             length(valid_idx) from fit_obj
##
## OUTPUT:
##   A list with:
##     $frac_state :
##        list with:
##          $pi1       : asymptotic mean fraction in state
##          $obs       : observed fraction in state
##          $var_CLT   : asymptotic variance of sample mean (≈ Var(π1_hat))
##          $se_CLT    : sqrt(var_CLT)
##          $z         : (obs - pi1) / se_CLT
##          $p_right   : P(Z >= z) under N(0,1) (obs > expected)
##     $transitions :
##        same Markov parameters used (pi1, pi11, pi01, pi00)
##
## NOTE:
##   - This function currently provides a CLT-based p-value only for
##     the fraction of time in the low-correlation state. Asymptotic
##     expectations for entry rate and run length are already given by
##     fit_obj$expected_asymptotic.

analytic_local_corr <- function(fit_obj,
                                T_eff = NULL) {
  
  print("Starting analytic Markov CLT inference (frac_state only)...")
  
  observed  <- fit_obj$observed
  expectedA <- fit_obj$expected_asymptotic
  
  ## Asymptotic quantities:
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
  ## Var(I_t) = pi1 (1 - pi1)
  ## theta = pi11 + pi00 - 1 (lag-1 autocorrelation of I_t)
  ## Var(π1_hat) ≈ Var(I_t) * (1 + theta) / (1 - theta) / T_eff
  var_CLT <- NA_real_
  se_CLT  <- NA_real_
  z       <- NA_real_
  p_right <- NA_real_
  
  if (is.finite(pi1_th) && pi1_th > 0 && pi1_th < 1 &&
      is.finite(pi11_th) && is.finite(pi00_th) &&
      T_eff > 0) {
    
    var_I   <- pi1_th * (1 - pi1_th)
    theta   <- pi11_th + pi00_th - 1  # corr(I_t, I_{t+1})
    
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
  
  print("Analytic Markov CLT inference complete.")
  
  list(
    frac_state = frac_state,
    transitions = transitions
  )
}