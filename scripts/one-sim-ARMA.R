library(tidyverse)
library(fda)
library(stats)
library(roll)

# generate one trend
n <- 100 # length of series to simulate
nb <- 25 # must be odd
s1 <- 1 # number of low-frequency components
s2 <- 3 # number of high-frequency components

# fourier basis functions
fb <- create.fourier.basis(rangeval = c(0, n), nbasis = nb, period = n)

# specify which frequencies will have nonzero coefficients
freq.sel <- 2*sample(5:((nb - 1)/2), s2, replace = F) - 1
freq.ix <- c(freq.sel, freq.sel + 1) |> sort()

# generate the coefficients
fb.coef <- rep(0, nb - 1)
fb.coef[1:(2*s1)] <- rnorm(2*s1, sd = 2)
fb.coef[freq.ix] <- rnorm(2*s2, sd = 0.25)

# construct the trend
x.state <- (eval.basis(1:n, basisobj = fb)[, -1] %*% fb.coef) |> scale()
plot(x.state, type = 'l')



# simulate noise series
ar.cor <- 0.5 # noise autocorrelation
# fixed sigma
sigma <- 0.7

# ARMA model coefficients

# Short Memory Coefs MA1
phi1 <- 0
phi2 <- 0
theta1 <- 0.25
theta2 <- 0

# # Slow decay Coefs with ARMA 21
# phi1 <- 0.8
# phi2 <- 0.1
# theta1 <- 0.3
# theta2 <- 0

ar.coefs <- c(phi1, phi2)   # AR coefs
ma.coefs <- c(theta1, theta2)         # MA coefs

noise <- arima.sim(
  model = list(order = c(2, 0, 2), ar = ar.coefs, ma = ma.coefs),
  n = n,
  rand.gen = function(t){ rnorm(t, mean = 0, sd = sigma) },
  n.start = 100
)

# add noise to drift
y <- x.state + noise

# Create lambda = observed variance / model implied variance

# fix S2_y to make it mean of sd of x.state over a rolling time k
# rolling average of x.state 
window <- 30 

loc_var <- roll_var(x.state, width = window, center = TRUE, min_obs = 2
)
s2_y_obs <- mean(loc_var, na.rm = TRUE)

# # window of lagmax
lag.max <- 50

# ARMAacf can return theoretical autocorrs for ARMA given noise variance
rho <- ARMAacf(ar = ar.coefs, ma = ma.coefs, lag.max = lag.max)
# converts the ARMA model into its infinite moving-average representation 
# so we can compute the theoretical stationary variance
# giving MA cofficients up to ag max - rerepre
psi <- ARMAtoMA(ar = ar.coefs, ma = ma.coefs, lag.max = lag.max)

# adds constant as the first variable to mak
psi_full <- c(1, psi)

# stationary variance gamma(0) for the ARMA process with noise sd = sigma
gamma0 <- (sigma^2) * sum(psi_full^2)

# convert rho(k) to autocovariance gamma(k)
gamma <- gamma0 * rho
# making k vector to iterate over
k <- 1:(window - 1)
# Equation 
sigma2_window_mean <- (1 / window^2) * (window * gamma[1] + 2 * sum((window - k) * gamma[k + 1]))  

sigma2_window_mean

plot(x.state, type = 'l')
lines(y, type = 'l')

lambda <- s2_y_obs / sigma2_window_mean

lambda

# showing decay structure of different coef groups
# candidates <- list(
#   short_MA1 = list(ar = numeric(0), ma = c(0.25)),
#   slow_ARMA21 = list(ar = c(0.80, 0.10), ma = c(0.30))
# )
# par(mfrow = c(2,2))
# for(nm in names(candidates)) {
#   ar <- candidates[[nm]]$ar
#   ma <- candidates[[nm]]$ma
#   rho <- ARMAacf(ar = ar, ma = ma, lag.max = lag.max)
#   plot(0:lag.max, rho, type="h", main = nm, xlab="lag", ylab="ACF")
#   abline(h=0)
# }
# par(mfrow = c(1,1))



