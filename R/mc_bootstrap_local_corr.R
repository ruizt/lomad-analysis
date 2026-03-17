#' Markov-chain bootstrap on the state process
#'
#' Parametric bootstrap on the binary state process \eqn{I_t} only, using a
#' 2-state Markov chain approximation derived from the asymptotic expectations
#' in [lomad_fit()]. More efficient than [lomad_test()] because
#' it does not re-simulate the full correlation process.
#'
#' @param fit_obj List returned by [lomad_fit()].
#' @param B Integer. Number of Markov-chain bootstrap replicates (default 1000).
#' @param T_sim Integer or NULL. Length of simulated chains. If `NULL`, uses
#'   `length(fit_obj$valid_idx)`.
#' @param boot_seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return A list containing:
#'   \describe{
#'     \item{expected_MC}{Bootstrap means under the Markov null: `entry_rate`,
#'       `mean_run_length`, `frac_state`.}
#'     \item{p_values_MC}{One-sided bootstrap p-values (obs >= null).}
#'     \item{transitions}{Markov parameters used: `pi1`, `pi0`, `pi11`,
#'       `pi01`, `pi00`.}
#'     \item{meta}{List with `B`, `T_sim`, and `boot_seed`.}
#'   }
#'
#' @export
mc_bootstrap_local_corr <- function(fit_obj,
                                    B = 1000,
                                    T_sim = NULL,
                                    boot_seed = NULL) {

  message("Starting Markov-chain bootstrap on I_t...")

  if (!is.null(boot_seed)) {
    set.seed(boot_seed)
  }

  observed  <- fit_obj$observed
  expectedA <- fit_obj$expected_asymptotic
  inputs    <- fit_obj$inputs

  pi1_th  <- expectedA$frac_state
  L_th    <- expectedA$mean_run_length
  pi11_th <- if (!is.na(L_th) && L_th > 1) 1 - 1 / L_th else NA_real_
  pi0_th  <- if (!is.na(pi1_th)) 1 - pi1_th else NA_real_

  ## pi01 = P(I_{t+1}=1 | I_t=0) implied by stationarity:
  ## pi1(1 - pi11) = pi0 * pi01  =>  pi01 = pi1(1-pi11) / (1-pi1)
  pi01_th <- if (!is.na(pi1_th) && !is.na(pi11_th) &&
                 pi1_th > 0 && pi1_th < 1) {
    pi1_th * (1 - pi11_th) / (1 - pi1_th)
  } else NA_real_

  pi00_th <- if (!is.na(pi01_th)) 1 - pi01_th else NA_real_

  if (is.null(T_sim)) {
    T_sim <- length(fit_obj$valid_idx)
  }

  if (any(!is.finite(c(pi1_th, pi11_th, pi01_th, pi00_th)))) {
    stop("Markov parameters (pi1, pi11, pi01, pi00) are not all finite.")
  }
  if (pi1_th <= 0 || pi1_th >= 1 ||
      pi11_th < 0 || pi11_th > 1 ||
      pi01_th < 0 || pi01_th > 1) {
    stop("Markov parameters fall outside [0,1]. Check fit/expected_asymptotic.")
  }

  ## Simulate one 2-state Markov chain of length T_sim
  simulate_I_chain <- function(T_sim, pi1, pi01, pi11) {
    I <- integer(T_sim)
    I[1] <- stats::rbinom(1, 1, pi1)
    if (T_sim > 1L) {
      for (t in 2:T_sim) {
        if (I[t - 1] == 1L) {
          I[t] <- stats::rbinom(1, 1, pi11)
        } else {
          I[t] <- stats::rbinom(1, 1, pi01)
        }
      }
    }
    I
  }

  ## Compute entry rate, mean run length, and frac_state from binary chain
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

  message(sprintf("Running Markov-chain bootstrap with %d replicates...", B))
  for (b in seq_len(B)) {
    I_sim <- simulate_I_chain(T_sim, pi1_th, pi01_th, pi11_th)
    stats_b <- compute_stats_from_I(I_sim)
    entry_boot[b]  <- stats_b$entry_rate
    runlen_boot[b] <- stats_b$mean_run_length
    frac_boot[b]   <- stats_b$frac_state
  }
  message("Markov-chain bootstrap complete.")

  entry_boot  <- entry_boot[is.finite(entry_boot)]
  runlen_boot <- runlen_boot[is.finite(runlen_boot)]
  frac_boot   <- frac_boot[is.finite(frac_boot)]

  expected_entry_MC <- if (length(entry_boot) > 0L) mean(entry_boot) else NA_real_
  expected_run_MC   <- if (length(runlen_boot) > 0L) mean(runlen_boot) else NA_real_
  expected_frac_MC  <- if (length(frac_boot) > 0L) mean(frac_boot) else NA_real_

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
    B         = B,
    T_sim     = T_sim,
    boot_seed = boot_seed
  )

  message("Markov-chain bootstrap finished.")

  list(
    expected_MC = expected_MC,
    p_values_MC = p_values_MC,
    transitions = transitions,
    meta        = meta
  )
}
