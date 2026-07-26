library(tidyverse)
library(fda)
library(stats)
library(roll)


# create first set of fourier coefficients
generate_fourier_coef <- function(nb, sd0 = 2, p=2.5) {
  K <- (nb - 1) / 2
  k <- rep(1:K, each = 2) #coefficient index vector
  sd_k <- sd0 / (k^p) 
  rnorm(2 * K, mean = 0, sd = sd_k)
}

generate_coef_pair <- function(nb, sd0 = 2, d = 1, p = 2.5, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  coef1 <- generate_fourier_coef(nb, sd0, p)
  
  # uniform sphere direction
  z <- rnorm(length(coef1))
  u <- z / sqrt(sum(z^2))  # unit length
  
  # ellipsoid axes (decay structure that is a function of distance)
  axes <- 1 / (1:length(coef1))^(p-0.1*d)
  
  # scale to ellipse
  dir <- u * axes
  
  # renormalize
  dir <- dir / sqrt(sum(dir^2))
  
  # shift by d
  coef2 <- coef1 + d * dir

  list(coef1 = coef1, coef2 = coef2)
}


# simulate latent trends
nb <- 25
n <- 500
d <- 1

coefs <- generate_coef_pair(nb = nb, d = d, seed=NULL)

coef1 <- coefs$coef1
coef2 <- coefs$coef2

fb <- create.fourier.basis(rangeval = c(0, n), nbasis = nb, period = n)

x1 <- as.numeric(eval.basis(1:n, fb)[, -1] %*% coef1) |> scale()
x2 <- as.numeric(eval.basis(1:n, fb)[, -1] %*% coef2) |> scale()

# plot(x1, type = 'l')
# plot(x2, type = 'l')


# ARMA noise parameters
# Low SNR 
lambda_target <- 0.5

# Slow decay Coefs with ARMA 21
ar.coefs <- c(0.8, 0.1)
ma.coefs <- c(0.3)

window <- 30 


# compute ARMA autocovariance structure
# ARMAacf can return theoretical autocorrs for ARMA given noise variance
rho <- ARMAacf(ar = ar.coefs, ma = ma.coefs, lag.max = (window - 1))
psi <- ARMAtoMA(ar = ar.coefs, ma = ma.coefs, lag.max = 200)
# contains all MA coeffs
psi_full <- c(1, psi)

# gamma0 when innovation variance is assumed 1,(because variance is 
# scaled linearly and we can separate arma shape from noise magnitude)
gamma0_unit <- sum(psi_full^2)          # gamma(0) when sigma^2 = 1
gamma_unit  <- gamma0_unit * rho        # gamma(k) when sigma^2 = 1
# making k vector to iterate over
k <- 1:(window - 1)
# Equation for variance of MA of length W of a stationary process
# plugging in gamma_unit which assumes innovation variance = 1
sigma2_window_mean_unit <- (1 / window^2) * (window * gamma_unit[1] + 2 * 
                                               sum((window - k) * gamma_unit[k + 1]))  

# fix S2_y to make it mean of sd of each series over a rolling time k
# rolling average of each series

loc_var_x1 <- roll_var(x1, width = window, center = TRUE, min_obs = 2
)
s2_y_obs_x1 <- mean(loc_var_x1, na.rm = TRUE)

loc_var_x2 <- roll_var(x2, width = window, center = TRUE, min_obs = 2
)
s2_y_obs_x2 <- mean(loc_var_x2, na.rm = TRUE)



#scale innovation sd for each series
sigma2_x1 <- s2_y_obs_x1 / (lambda_target * sigma2_window_mean_unit)
sigma_x1  <- sqrt(sigma2_x1)

sigma2_x2 <- s2_y_obs_x2 / (lambda_target * sigma2_window_mean_unit)
sigma_x2  <- sqrt(sigma2_x2)

noise1 <- arima.sim(
  model = list(order = c(length(ar.coefs), 0, length(ma.coefs)), ar = ar.coefs, ma = ma.coefs),
  n = n,
  rand.gen = function(t){ rnorm(t, mean = 0, sd = sigma_x1) },
  n.start = 100)

noise2 <- arima.sim(
  model = list(order = c(length(ar.coefs), 0, length(ma.coefs)), ar = ar.coefs, ma = ma.coefs),
  n = n,
  rand.gen = function(t){ rnorm(t, mean = 0, sd = sigma_x2) },
  n.start = 100)


# add noise to drift
y1 <- x1 + noise1
y2 <- x2 + noise2

plot(y1, type = 'l', col = "black", ylim=range(c(y1,y2)))
lines(y2, type = 'l', col = "blue")


