library(tidyverse)
library(fda)
library(ggplot2)

# create first set of fourier coefficients
generate_fourier_coef <- function(nb, sd0 = 2, p=2.5) {
  K <- (nb - 1) / 2
  k <- rep(1:K, each = 2) #coefficient index vector
  sd_k <- sd0 / (k^p) 
  rnorm(2 * K, mean = 0, sd = sd_k)
}


#### UNIT SPHERE METHOD ####
# 
# # generate pair of coefficients
# generate_coef_pair <- function(nb, sd0 = 2, d = 1, p=2.5, seed = NULL) {
#   if (!is.null(seed)) set.seed(seed)
# 
#   # first coefficient vector
#   coef1 <- generate_fourier_coef(nb, sd0)
# 
#   # second coefficient vector at distance d
#   z <- rnorm(length(coef1)) # random direction on unit sphere
#   direction <- z / sqrt(sum(z^2)) #direction u vect
#   coef2 <- coef1 + d * direction
# 
#   list(coef1 = coef1, coef2 = coef2)
# }


#### SPHERE TO ELLIPSE TRANSFORM METHOD ####
# slightly different decay structures each time

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



#### Test Euclidean Distance ####
coefs <- generate_coef_pair(25, d=2)
sqrt(sum((coefs$coef2 - coefs$coef1)^2))




#### Plot Series Against Each Other ####
nb <- 25
n <- 500
d <- 2

coefs <- generate_coef_pair(nb = nb, d = d, seed=NULL)

coef1 <- coefs$coef1
coef2 <- coefs$coef2

fb <- create.fourier.basis(rangeval = c(0, n), nbasis = nb, period = n)

Phi <- eval.basis(1:n, fb)[, -1]

x1 <- as.numeric(Phi %*% coef1)
x2 <- as.numeric(Phi %*% coef2)


plot( x1, type = "l", col = "black", ylim=range(c(x1,x2)))
lines(x2, col = "blue")
legend( "topright", legend = c("Coef 1", "Coef 2"), col = c("black", "blue"), lwd=2,bty="n")

