library(tidyverse)

# load lomad package in place
library(lomad)

# read in data
simdata <- read_rds('example-data/sim-example.rds')
str(simdata)

# compute moving averages and rolling correlation
df <- simdata$`1`
lfit <- lomad_fit(df$x1, df$x2, q = 10, h = 30)

# moving averages
x1.ma <- lfit$ma1
x2.ma <- lfit$ma2
plot(x1.ma, type = 'l', col = 'red')
lines(x2.ma, col = 'blue')

# true trends (simulated data only)
plot(df$x1_trend, type = 'l', col = 'red', lty = 1)
lines(df$x2_trend, col = 'blue', lty = 2)

# correlation series
rt <- lfit$R
plot(rt, type = 'l')
