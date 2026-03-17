devtools::load_all()

# generate trends
trends <- make_trends_dist(n = 500, d = 5, seed = 32026)
plot(trends$x1, type = "l", col = 'red',
     xlab = "t", ylab = expression(mu[t]))
lines(trends$x2, col = "blue")

# add ARMA noise
sim <- add_noise(trends,
                 h             = 30,
                 lambda_target = 4,
                 scale         = 5,
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

# local moving average decoupling
out <- lomad(sim$y1, sim$y2, q = 30, h = 50, B = 1000, seed = 31726)
out$p_values
out$observed
out$expected

# plot
plot_lomad_fit(sim$y1, sim$y2, out)

# plot smoothed trends only
plot(out$ma1, type = "l", col = "red")
lines(out$ma2, type = "l", col = "blue")
lines(out$trend_hat, type = 'l', col = 'grey')
