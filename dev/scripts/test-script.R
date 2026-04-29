devtools::load_all()

# generate trends
trends <- sim_trends(n = 500, d = 2, method = "dist", seed = 32026)
plot(trends$x1, type = "l", col = 'red',
     xlab = "t", ylab = expression(mu[t]))
lines(trends$x2, col = "blue")

# add ARMA noise
sim <- sim_noise_pair(trends,
                      h             = 30,
                      lambda_target = 4,
                      scale         = 1,
                      order         = c(2, 1),
                      s             = 100,
                      n_start       = 30,
                      seed          = 31726)

# summary of added noise processes
sim$noise

# plot observed series with true trends overlaid
plot(sim$y1, type = "l", col = 'red',
     xlab = 't', ylab = expression(X[t]))
lines(sim$y2, col = "blue")
lines(sim$x1, col = "darkgrey", lwd = 2)
lines(sim$x2, col = "darkgrey", lwd = 2)

# local moving average decoupling — state method with bootstrap test
out <- lomad(sim$y1, sim$y2, method = "state",
             test_method = "boot",
             q = 10, h = 50, B = 100,
             seed = 31726, ncores = parallel::detectCores() - 1,
             verbose = TRUE)
out$test$p_values

# plot
lomad_plot(out$fit, x1 = sim$y1, x2 = sim$y2)

# plot smoothed trends only
plot(out$fit$ma1, type = "l", col = "red")
lines(out$fit$ma2, type = "l", col = "blue")
lines(out$fit$trend_hat, type = 'l', col = 'grey')
