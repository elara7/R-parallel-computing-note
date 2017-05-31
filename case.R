library(risksetROC)
library(MASS)
library(dynpred)
library(mvtnorm)
library(Kendall)
library(smooth)
library(data.table)

mu1 <- mu2 <- 0
sd1 <- sd2 <- 1
rho <- -0.7
mu <- c(mu1,mu2)
sigma <- matrix(c(sd1^2, sd1*sd2*rho, sd1*sd2*rho, sd2^2), 2)

FPDynamic <- function(data, rho) {
  c <- data[1]
  logt <- data[2]
  pmvnorm(c(-Inf, -Inf), c(-c, -logt),mean = mu, 
          sigma = matrix(c(1, rho, rho, 1), 2)) / pnorm(- logt)
}

TPIncident <- function(data, rho) {
  c <- data[1]
  logt <- data[2]
  pnorm( (rho * logt - c) / sqrt(1 - rho^2) )
}

AUCNormTrue <- function(logt, rho) {
  # Calculate the true AUC value under certain distribution/
  # Bi-Norm(0,0,1,1,rho)
  c <- seq(-5, 5, 0.01)
  tp <- pnorm( (rho * logt - c) / sqrt(1 - rho^2) )
  fp <- apply(cbind(c, rep(logt, length(c))), 1, FPDynamic, rho = rho)
  # true value of AUC
  dFP <- abs(fp[-1] - fp[-length(fp)])
  aTP <- 0.5 * (tp[-1] + tp[-length(tp)])
  sum(dFP * aTP)
}

seq.ltime <- seq(-2, 2, 0.5)
seq.time <- exp(seq.ltime)

theoAUC <- sapply(seq.ltime, AUCNormTrue, rho = -0.7)


# choose the distb for 20% censoring
pnorm(0, mean = 1.1901, sd = 1.414)
mean.20 <- 1.1901
# 40% censoring
pnorm(0, mean = 0.3584, sd = 1.414)
mean.40 <- 0.3584

set.seed(2017)
N <- 200*1000
SIM1 <- data.frame(mvrnorm(N, mu = mu, Sigma = sigma))
SIM1$M <- SIM1$X1
SIM1$X1 <- NULL
SIM1$censored.20 <- rnorm(N, mean.20, 1)
SIM1$censored.40 <- rnorm(N, mean.40, 1)
SIM1$status.20 <- SIM1$censored.20 > SIM1$X2
SIM1$status.40 <- SIM1$censored.40 > SIM1$X2
SIM1$T.20 <- exp(apply(cbind(SIM1$censored.20, SIM1$X2), 1, min))
SIM1$T.40 <- exp(apply(cbind(SIM1$censored.40, SIM1$X2), 1, min))
SIM1$IND <- rep(1:1000, each = 200)
SIM <- split(SIM1, SIM1$IND)
