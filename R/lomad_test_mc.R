#' Markov-chain bootstrap test for local correlation decoupling
#'
#' Parametric bootstrap on the binary state process \eqn{I_t} only, using a
#' 2-state Markov chain approximation derived from [lomad_fit()]. More
#' efficient than [lomad_test_boot()] because it does not re-simulate the full
#' correlation process.
#'
#' @param fit List returned by [lomad_fit()].
#' @param B Integer. Number of bootstrap replicates (default 1000).
#' @param T_sim Integer or NULL. Length of simulated chains. If `NULL`,
#'   defaults to `length(fit$valid_idx)`.
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return A named list with the same structure as [lomad_test_boot()]:
#'   \describe{
#'     \item{null_model}{Null model from `fit`, augmented with Markov
#'       transition parameters (`transitions`).}
#'     \item{observed}{Observed statistics from `fit`.}
#'     \item{expected}{Bootstrap means under the Markov null.}
#'     \item{expected_asymptotic}{Asymptotic expectations from `fit`.}
#'     \item{p_values}{One-sided bootstrap p-values (P(stat >= observed)).}
#'     \item{inputs}{`n`, `max_pq`, `B`, `T_sim`, `seed`.}
#'     \item{trend_hat, ma1, ma2, R, thresholds, I, valid_idx}{Passed through
#'       from `fit` for plotting and downstream use.}
#'   }
#'
#' @export
lomad_test_mc <- function(fit,
                           B     = 1000,
                           T_sim = NULL,
                           seed  = NULL) {

  if (!is.null(seed)) set.seed(seed)

  observed  <- fit$observed
  expectedA <- fit$expected_asymptotic

  # --- Markov transition parameters from asymptotic expectations ---

  pi1_th  <- expectedA$frac_state
  L_th    <- expectedA$mean_run_length
  pi11_th <- if (!is.na(L_th) && L_th > 1) 1 - 1 / L_th else NA_real_
  pi0_th  <- if (!is.na(pi1_th)) 1 - pi1_th else NA_real_
  pi01_th <- if (!is.na(pi1_th) && !is.na(pi11_th) &&
                 pi1_th > 0 && pi1_th < 1)
    pi1_th * (1 - pi11_th) / (1 - pi1_th) else NA_real_
  pi00_th <- if (!is.na(pi01_th)) 1 - pi01_th else NA_real_

  if (any(!is.finite(c(pi1_th, pi11_th, pi01_th, pi00_th))))
    stop("Markov transition parameters are not all finite. Check fit$expected_asymptotic.")
  if (pi1_th <= 0 || pi1_th >= 1 ||
      pi11_th < 0 || pi11_th > 1 ||
      pi01_th < 0 || pi01_th > 1)
    stop("Markov transition parameters fall outside [0, 1].")

  if (is.null(T_sim)) T_sim <- length(fit$valid_idx)

  # --- Helper: simulate one 2-state Markov chain ---

  sim_chain <- function() {
    I    <- integer(T_sim)
    I[1] <- stats::rbinom(1, 1, pi1_th)
    for (t in 2:T_sim)
      I[t] <- stats::rbinom(1, 1, if (I[t - 1] == 1L) pi11_th else pi01_th)
    I
  }

  # --- Helper: state statistics from one simulated chain ---

  chain_stats <- function(I_v) {
    m          <- length(I_v)
    frac_state <- mean(I_v)

    steps01    <- sum(I_v[-m] == 0L)
    entry_rate <- if (steps01 > 0L)
      sum(I_v[-1L] == 1L & I_v[-m] == 0L) / steps01 else NA_real_

    runs            <- rle(I_v)
    run_lengths     <- runs$lengths[runs$values == 1L]
    n_entries_sim   <- length(run_lengths)
    mean_run_length <- if (n_entries_sim > 0L) mean(run_lengths) else NA_real_

    list(entry_rate      = entry_rate,
         mean_run_length = mean_run_length,
         frac_state      = frac_state,
         n_entries       = as.numeric(n_entries_sim))
  }

  # --- Bootstrap loop ---

  boot <- lapply(seq_len(B), function(b) chain_stats(sim_chain()))

  entry_v <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "entry_rate"))
  run_v   <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "mean_run_length"))
  frac_v  <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "frac_state"))
  nent_v  <- Filter(is.finite, vapply(boot, `[[`, numeric(1), "n_entries"))

  # --- Summaries ---

  expected <- list(
    entry_rate      = if (length(entry_v)) mean(entry_v) else NA_real_,
    mean_run_length = if (length(run_v))   mean(run_v)   else NA_real_,
    frac_state      = if (length(frac_v))  mean(frac_v)  else NA_real_,
    n_entries       = if (length(nent_v))  mean(nent_v)  else NA_real_
  )

  p_values <- list(
    entry_rate      = if (!is.na(observed$entry_rate)      && length(entry_v))
      mean(entry_v >= observed$entry_rate)      else NA_real_,
    mean_run_length = if (!is.na(observed$mean_run_length) && length(run_v))
      mean(run_v   >= observed$mean_run_length) else NA_real_,
    frac_state      = if (!is.na(observed$frac_state)      && length(frac_v))
      mean(frac_v  >= observed$frac_state)      else NA_real_,
    n_entries       = if (!is.na(observed$n_entries)       && length(nent_v))
      mean(nent_v  >= observed$n_entries)       else NA_real_
  )

  message(sprintf(
    "B = %d | p: entry_rate = %s, run_length = %s, frac_state = %s, n_entries = %s",
    B,
    signif(p_values$entry_rate, 3),
    signif(p_values$mean_run_length, 3),
    signif(p_values$frac_state, 3),
    signif(p_values$n_entries, 3)
  ))

  null_model             <- fit$null_model
  null_model$transitions <- list(pi1  = pi1_th,  pi0  = pi0_th,
                                  pi11 = pi11_th, pi01 = pi01_th, pi00 = pi00_th)

  list(
    null_model          = null_model,
    observed            = observed,
    expected            = expected,
    expected_asymptotic = expectedA,
    p_values            = p_values,
    inputs              = list(n = fit$inputs$n, max_pq = fit$inputs$max_pq,
                               B = B, T_sim = T_sim, seed = seed),
    trend_hat  = fit$trend_hat,
    ma1        = fit$ma1,
    ma2        = fit$ma2,
    R          = fit$R,
    thresholds = fit$thresholds,
    I          = fit$I,
    valid_idx  = fit$valid_idx
  )
}
