library(lomad)

#### old generate_lambda ####
# "decouple_rate" as input

generate_lambda_gauss <- function(n,
                            decouple_rate = 0.1,  
                            decouple_strength = 1, # intesity of decoupling spike (0 to 1)
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n_events <- floor(1 / decouple_rate) #number of decouplings in series
  gap <- n / n_events # distance between centers of each decoupling
  decouple_length <- round(gap * decouple_rate)  #length of each decoupling; wider as rate increases
  
  #space out decoupling events
  #gives midpoints of each decoupling event
  events <- round(seq(from = gap / 2, by = gap, length.out = n_events))
  
  t <- 1:n
  lambda <- rep(1, n) #initialize weights vector
  
  # create decouplings
  for (i in events) {
    bump <- dnorm(t, mean = i, sd = decouple_length / 4) # normal curve centered at i
    bump <- bump / max(bump) # normalize so peak is 1
    
    # create dips in lambda by subtracting bump
    # made function of rate and strength; lower values of either => shallower dips
    lambda <- lambda - decouple_rate * decouple_strength * bump 
  }
  
  pmax(0, pmin(1, lambda)) #ensure lambda within [0,1]
}




#### generate gamma decouplings ####

# create function for gamma decoupling event
# want gamma bumps centered at event midpoints
# make scale a function of decoupling length ?

#also change to make n_events input instead of "decouple_rate"

generate_lambda_gamma <- function(n,
                                  n_events = 5,  # number of decoupling events
                                  decouple_length = NULL, # length of each decoupling; NULL => compute from n_events
                                  decouple_strength = 1, # intesity of decoupling spike
                                  shape = 2, # shape parameter of gamma
                                  sigma_scale = 0.25,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  delta <- 1 / n_events # fraction of timeline each decoupling event owns (previously ~deouple_rate)
  gap <- n / n_events # distance between centers of each decoupling
 
 # use decouple_length input if given, otherwise derive from n_events
  if (is.null(decouple_length)){
    decouple_length <- round(gap * delta)  #length of each decoupling; narrower as n_events increases
  }
  
  #space out decoupling events
  #gives midpoints of each decoupling event
  events <- round(seq(from = gap / 2, by = gap, length.out = n_events))
  
  t <- 1:n
  lambda <- rep(1, n) #initialize weights vector
  
  # create decouplings using a gamma curve
  for (i in events) {
    
    bump <- dgamma(t-i, shape = shape, scale = decouple_length * sigma_scale) # normal curve centered at i
    bump <- bump / max(bump) # normalize so peak is 1
    
    # create dips in lambda by subtracting bump
    # made function of delta and strength; lower values of either => shallower dips
    lambda <- lambda - delta * decouple_strength * bump 
  }
  
  pmax(0, pmin(1, lambda)) #ensure lambda within [0,1]
}



#### PLOT SERIES ####
# create series
n <- 500
nb <- 25
decouple_rate <- 0.5
n_events <- 5
decouple_length = NULL
sigma_scale = 0.5
decouple_strength <- 3
shape <- 2
seed <- NULL

## fixed distance coefficents ##
coefs <- generate_coef_pair(nb,d=1)
coef1 <- coefs$coef1
coef2 <- coefs$coef2

lambda <- generate_lambda_gamma(n, n_events, decouple_length, decouple_strength, shape, sigma_scale)

fb  <- fda::create.fourier.basis(rangeval = c(0, n), nbasis = nb, period = n)
Phi <- fda::eval.basis(1:n, fb)[, -1]

x1 <- as.numeric(Phi %*% coef1)
x2 <- as.numeric(Phi %*% coef2)
x_mixed <- lambda*x1 + (1-lambda)*x2

plot( x1, type = "l", col = "black", ylim=range(c(x1,x2)))
lines (x2, col="gray")
lines(x_mixed, col = "blue")
