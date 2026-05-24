# Internal utilities for estimate_*() functions
#
# .variogram_ar1()    — AR(1) estimation from lag-1/2 variogram
# .variogram_ar()     — general AR(p) from multi-lag variogram + Yule-Walker
# .yule_walker()      — solve Yule-Walker equations
# .fit_arma()         — AIC-based ARMA fitting (retained for future use)
# .select_arma()      — AIC-based ARMA order selection (used by lomad_test_identity)
# .smooth_noise_var() — variance of MA-smoothed noise process (used by lomad_test_identity)


# AR(1) estimation from the lag-1 and lag-2 variogram.
.variogram_ar1 <- function(resid) {
  n  <- length(resid)
  d1 <- diff(resid)
  d2 <- resid[3:n] - resid[1:(n - 2L)]

  V1 <- mean(d1^2) / 2
  V2 <- mean(d2^2) / 2

  phi <- V2 / V1 - 1
  phi <- min(max(phi, 0.01), 0.99)

  list(
    ar     = phi,
    ma     = numeric(0),
    sigma2 = V2,
    order  = c(1L, 0L, 0L)
  )
}


# General AR(p) estimation from multi-lag variogram + Yule-Walker.
.variogram_ar <- function(resid, p = NULL, p_max = 5L, L = 20L) {
  n <- length(resid)
  if (L >= n - 1L) L <- n - 2L
  if (is.null(p)) {
    max_order <- p_max
  } else {
    max_order <- p
  }
  if (L < max_order + 1L) L <- max_order + 1L

  V <- numeric(L)
  for (l in seq_len(L)) {
    dl <- resid[(l + 1L):n] - resid[1:(n - l)]
    V[l] <- mean(dl^2) / 2
  }

  plateau_start <- max(1L, ceiling(L / 2))
  gamma0 <- mean(V[plateau_start:L])
  gamma_hat <- gamma0 - V

  if (!is.null(p)) {
    yw <- .yule_walker(gamma0, gamma_hat, p)
    return(list(
      ar     = yw$ar,
      ma     = numeric(0),
      sigma2 = yw$sigma2,
      order  = c(p, 0L, 0L)
    ))
  }

  best_bic <- Inf
  best_p   <- 1L
  best_yw  <- NULL

  for (pp in seq_len(p_max)) {
    yw  <- .yule_walker(gamma0, gamma_hat, pp)
    if (yw$sigma2 <= 0) next
    bic <- n * log(yw$sigma2) + pp * log(n)
    if (bic < best_bic) {
      best_bic <- bic
      best_p   <- pp
      best_yw  <- yw
    }
  }

  if (is.null(best_yw)) {
    best_yw <- .yule_walker(gamma0, gamma_hat, 1L)
    best_p  <- 1L
  }

  list(
    ar     = best_yw$ar,
    ma     = numeric(0),
    sigma2 = best_yw$sigma2,
    order  = c(best_p, 0L, 0L)
  )
}


# Solve Yule-Walker equations for AR(p).
.yule_walker <- function(gamma0, gamma_hat, p) {
  acvf <- c(gamma0, gamma_hat[seq_len(p)])
  Gamma_mat <- stats::toeplitz(acvf[seq_len(p)])
  gamma_vec <- gamma_hat[seq_len(p)]

  phi <- tryCatch(
    solve(Gamma_mat, gamma_vec),
    error = function(e) rep(0, p)
  )

  sigma2 <- gamma0 - sum(phi * gamma_vec)
  sigma2 <- max(sigma2, 1e-10)

  list(ar = as.numeric(phi), sigma2 = sigma2)
}


# AIC-based ARMA order selection and fitting.
# Retained for potential future use; not called by the default pipeline.
.fit_arma <- function(resid, max_pq) {
  best_aic   <- Inf
  best_fit   <- NULL
  best_order <- c(0L, 0L, 0L)

  for (p in 0:max_pq) {
    for (q in 0:max_pq) {
      if (p == 0L && q == 0L) next
      converged <- TRUE
      fit <- withCallingHandlers(
        try(stats::arima(resid, order = c(p, 0L, q),
                         optim.control = list(maxit = 500L)),
            silent = TRUE),
        warning = function(w) {
          msg <- conditionMessage(w)
          if (grepl("convergence problem|NaN|log", msg))
            invokeRestart("muffleWarning")
          if (grepl("convergence problem", msg))
            converged <<- FALSE
        }
      )
      if (inherits(fit, "try-error") || !converged) next
      if (is.na(fit$aic) || !is.finite(fit$aic)) next
      if (is.null(fit$sigma2) || fit$sigma2 <= 0) next
      if (fit$aic < best_aic) {
        best_aic   <- fit$aic
        best_fit   <- fit
        best_order <- c(p, 0L, q)
      }
    }
  }

  if (is.null(best_fit)) {
    best_fit   <- stats::arima(resid, order = c(0L, 0L, 0L))
    best_order <- c(0L, 0L, 0L)
  }

  coef <- best_fit$coef
  list(
    ar     = as.numeric(coef[grepl("^ar", names(coef))]),
    ma     = as.numeric(coef[grepl("^ma", names(coef))]),
    sigma2 = best_fit$sigma2,
    order  = best_order,
    resid  = as.numeric(stats::residuals(best_fit))
  )
}


# [LEGACY] AIC-based ARMA order selection — used by the state-process
# pipeline (.lomad_fit_state, .lomad_fit_blocks) and lomad_test_identity.
.select_arma <- function(resid, max_pq) {
  best_aic   <- Inf
  best_fit   <- NULL
  best_order <- c(0L, 0L, 0L)

  for (p in 0:max_pq) {
    for (qi in 0:max_pq) {
      if (p == 0L && qi == 0L) next
      converged <- TRUE
      fit <- withCallingHandlers(
        try(stats::arima(resid, order = c(p, 0L, qi),
                         optim.control = list(maxit = 500L)), silent = TRUE),
        warning = function(w) {
          if (grepl("convergence problem", conditionMessage(w))) {
            converged <<- FALSE
            invokeRestart("muffleWarning")
          }
        }
      )
      if (inherits(fit, "try-error") || is.na(fit$aic) || !converged) next
      if (fit$aic < best_aic) {
        best_aic   <- fit$aic
        best_fit   <- fit
        best_order <- c(p, 0L, qi)
      }
    }
  }
  if (is.null(best_fit)) {
    best_fit   <- stats::arima(resid, order = c(0L, 0L, 0L))
    best_order <- c(0L, 0L, 0L)
  }
  coef <- best_fit$coef
  list(order  = best_order,
       ar     = as.numeric(coef[grepl("^ar", names(coef))]),
       ma     = as.numeric(coef[grepl("^ma", names(coef))]),
       sigma2 = best_fit$sigma2)
}


# Variance of the MA(q)-smoothed noise process.
# Used by lomad_test_identity().
.smooth_noise_var <- function(fit, q) {
  ar     <- fit$ar
  ma     <- fit$ma
  sigma2 <- fit$sigma2
  max_lag <- q - 1L

  if (length(ar) == 0L && length(ma) == 0L) {
    return(sigma2 / q)
  }

  psi    <- c(1, stats::ARMAtoMA(ar = ar, ma = ma,
                                  lag.max = max(200L, max_lag)))
  gamma0 <- sigma2 * sum(psi^2)
  acf    <- stats::ARMAacf(ar = ar, ma = ma, lag.max = max_lag)
  gamma  <- gamma0 * as.numeric(acf)

  weights <- q - seq_len(max_lag)
  (q * gamma[1L] + 2 * sum(weights * gamma[-1L])) / q^2
}
