## mc_bootstrap_local_corr.R
##
## Parametric bootstrap on the *state process* I_t only, using a
## 2-state Markov chain approximation derived from the asymptotic
## expectations in fit_null_model_local_corr().
##
## INPUTS:
##   fit_obj   : list returned by fit_null_model_local_corr()
##   B         : number of Markov-chain bootstrap replicates (default 1000)
##   T_sim     : length of simulated chains; if NULL, uses length(valid_idx)
##   boot_seed : (optional) integer RNG seed for reproducible bootstrap
##
## OUTPUT:
##   A list with:
##     $expected_MC : bootstrap means under the Markov null:
##                    $entry_rate, $mean_run_length, $frac_state
##     $p_values_MC : one-sided bootstrap p-values (obs >= null)
##     $transitions : list with Markov parameters used (pi1, pi11, pi01, pi00)
##     $meta        : list with B, T_sim, boot_seed

mc_bootstrap_local_corr <- function(fit_obj,
                                              B = 1000,
                                              T_sim = NULL,
                                              boot_seed = NULL) {
  
  print("Starting Markov-chain bootstrap on I_t...")
  
  if (!is.null(boot_seed)) {
    set.seed(boot_seed)
  }
  
  observed  <- fit_obj$observed
  expectedA <- fit_obj$expected_asymptotic
  inputs    <- fit_obj$inputs
  
  ## Asymptotic quantities:
  pi1_th  <- expectedA$frac_state       # π1
  L_th    <- expectedA$mean_run_length  # ≈ 1/(1-π11)
  ## We recompute π11, π01 from π1 and L_th to enforce Markov consistency
  pi11_th <- if (!is.na(L_th) && L_th > 1) 1 - 1 / L_th else NA_real_
  pi0_th  <- if (!is.na(pi1_th)) 1 - pi1_th else NA_real_
  
  ## π01 = P(I_{t+1}=1 | I_t=0) implied by stationarity:
  ## π1(1 - π11) = π0 π01  =>  π01 = π1(1-π11) / (1-π1)
  pi01_th <- if (!is.na(pi1_th) && !is.na(pi11_th) &&
                 pi1_th > 0 && pi1_th < 1) {
    pi1_th * (1 - pi11_th) / (1 - pi1_th)
  } else NA_real_
  
  pi00_th <- if (!is.na(pi01_th)) 1 - pi01_th else NA_real_
  
  if (is.null(T_sim)) {
    ## Use number of valid correlation windows as effective chain length
    T_sim <- length(fit_obj$valid_idx)
  }
  
  ## Sanity checks
  if (any(!is.finite(c(pi1_th, pi11_th, pi01_th, pi00_th)))) {
    stop("Markov parameters (pi1, pi11, pi01, pi00) are not all finite.")
  }
  if (pi1_th <= 0 || pi1_th >= 1 ||
      pi11_th < 0 || pi11_th > 1 ||
      pi01_th < 0 || pi01_th > 1) {
    stop("Markov parameters fall outside [0,1]. Check fit/expected_asymptotic.")
  }
  
  ## Helper: simulate one 2-state Markov chain of length T_sim
  simulate_I_chain <- function(T_sim, pi1, pi01, pi11) {
    I <- integer(T_sim)
    ## start from stationary distribution
    I[1] <- rbinom(1, 1, pi1)
    if (T_sim > 1L) {
      for (t in 2:T_sim) {
        if (I[t - 1] == 1L) {
          I[t] <- rbinom(1, 1, pi11)
        } else {
          I[t] <- rbinom(1, 1, pi01)
        }
      }
    }
    I
  }
  
  ## Helper: compute stats from a binary chain I_t
  compute_stats_from_I <- function(I_vec) {
    n <- length(I_vec)
    if (n <= 1L) {
      return(list(entry_rate = NA_real_,
                  mean_run_length = NA_real_,
                  frac_state = mean(I_vec)))
    }
    
    frac_state <- mean(I_vec)
    
    I_prev <- I_vec[-n]
    I_curr <- I_vec[-1L]
    
    steps01 <- sum(I_prev == 0L)
    entries <- sum(I_prev == 0L & I_curr == 1L)
    entry_rate <- if (steps01 > 0L) entries / steps01 else NA_real_
    
    run_lengths <- integer(0)
    if (any(I_vec == 1L)) {
      rle_I <- rle(I_vec)
      run_lengths <- rle_I$lengths[rle_I$values == 1L]
    }
    mean_run_length <- if (length(run_lengths) > 0L) mean(run_lengths) else NA_real_
    
    list(entry_rate = entry_rate,
         mean_run_length = mean_run_length,
         frac_state = frac_state)
  }
  
  entry_boot  <- numeric(B)
  runlen_boot <- numeric(B)
  frac_boot   <- numeric(B)
  
  cat("Running Markov-chain bootstrap with", B, "replicates:\n")
  for (b in seq_len(B)) {
    if (b == 1) {
      cat("  Iteration:")
    }
    if (b %% 10 == 1) {
      cat("\n   ")
    }
    cat(b, "")
    
    I_sim <- simulate_I_chain(T_sim, pi1_th, pi01_th, pi11_th)
    stats_b <- compute_stats_from_I(I_sim)
    entry_boot[b]  <- stats_b$entry_rate
    runlen_boot[b] <- stats_b$mean_run_length
    frac_boot[b]   <- stats_b$frac_state
  }
  cat("\nMarkov-chain bootstrap complete.\n")
  
  entry_boot  <- entry_boot[is.finite(entry_boot)]
  runlen_boot <- runlen_boot[is.finite(runlen_boot)]
  frac_boot   <- frac_boot[is.finite(frac_boot)]
  
  expected_entry_MC <- if (length(entry_boot) > 0L) mean(entry_boot) else NA_real_
  expected_run_MC   <- if (length(runlen_boot) > 0L) mean(runlen_boot) else NA_real_
  expected_frac_MC  <- if (length(frac_boot) > 0L) mean(frac_boot) else NA_real_
  
  ## One-sided p-values: obs >= null
  p_entry_MC <- if (!is.na(observed$entry_rate) && length(entry_boot) > 0L) {
    mean(entry_boot >= observed$entry_rate)
  } else NA_real_
  
  p_run_MC <- if (!is.na(observed$mean_run_length) && length(runlen_boot) > 0L) {
    mean(runlen_boot >= observed$mean_run_length)
  } else NA_real_
  
  p_frac_MC <- if (!is.na(observed$frac_state) && length(frac_boot) > 0L) {
    mean(frac_boot >= observed$frac_state)
  } else NA_real_
  
  expected_MC <- list(
    entry_rate      = expected_entry_MC,
    mean_run_length = expected_run_MC,
    frac_state      = expected_frac_MC
  )
  
  p_values_MC <- list(
    entry_rate      = p_entry_MC,
    mean_run_length = p_run_MC,
    frac_state      = p_frac_MC
  )
  
  transitions <- list(
    pi1  = pi1_th,
    pi0  = pi0_th,
    pi11 = pi11_th,
    pi01 = pi01_th,
    pi00 = pi00_th
  )
  
  meta <- list(
    B        = B,
    T_sim    = T_sim,
    boot_seed = boot_seed
  )
  
  print("Markov-chain bootstrap finished.")
  
  list(
    expected_MC = expected_MC,
    p_values_MC = p_values_MC,
    transitions = transitions,
    meta = meta
  )
}