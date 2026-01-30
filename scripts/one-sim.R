library(tidyverse)
library(fda)

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
fb.coef <- c(rnorm(n = 4, mean = 1, sd = 0.5), rep(0, 20))

# construct the trend
x.state <- (eval.basis(1:n, basisobj = fb)[, -1] %*% fb.coef) |> scale()
plot(x.state, type = 'l')

# simulate noise series
ar.cor <- 0.5 # noise autocorrelation
sigma <- 0.5 # noise standard deviation
noise <- arima.sim(model = list(order = c(1, 0, 0), ar = ar.cor),
                      n = n,
                      rand.gen = function(t){rnorm(t, mean = 0, sd = sigma)},
                      n.start = 100)

# add noise to drift
y <- x.state + noise

plot(x.state, type = 'l')
lines(y, type = 'l')
