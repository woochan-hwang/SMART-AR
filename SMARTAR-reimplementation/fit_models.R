# Reimplementation and extension of SMART with adaptive randomisation [Biometrics, Cheung 2015]
# Woochan H. 2020, supervised by Dr David Robertson & Dr Sofia VIllar

# Load dependencies ---------------------
library(magrittr)


# Q-functions ---------------------------
# Saturated Q2 model used as gold standard to simulate patient outcome
Q2sat <- function(a1,resp,a2, theta_sat) {
  dim(theta_sat) <- c(8,1)
  x <- cbind(1, a1, a2, resp, a1*a2, a1*resp, a2*resp, a1*a2*resp)
  return (x %*% theta_sat)
}

# Q2 parameterised to allow PTW & PTW-m strategy
Q2 <- function(a1,resp,a2,b0,b1,b2,b3,g1,g2,g3) {
  return (b0 + b1*a1 + b2*a2 + b3*a1*a2 + g1*resp + g2*resp*(1-a1)*a2 + g3*resp*a1*(1-a2))
}

# pseudo-outcome used for Q1
Q2max <- function(a1, resp, b0, b1, b2, b3, g1, g2, g3) {
  Q20 <- Q2(a1, resp, 0, b0, b1, b2, b3, g1, g2, g3)
  Q21 <- Q2(a1, resp, 1, b0, b1, b2, b3, g1, g2, g3)
  if (Q20 >= Q21)  return( c(0,Q20) )
  else { return ( c(1,Q21) ) }
}

# not used in Qlearn because p = response rate under a1 is not available
Q1 <- function(a1, b0, b1, b2, b3, g1, g2, g3, p) {
  val <- (1-p) * Q2max(a1, 0, b0, b1, b2, b3, g1, g2, g3)[2] + p * Q2max(a1, 1, b0, b1, b2, b3, g1, g2, g3)[2]
  return (val)
}


# Fit Q-functions ------------------------------------------
FitQ2 <- function(foo) {
  z1 <- a1 <- foo[, 1]
  z2 <- a2 <- foo[, 3]
  z3 <- a1*a2
  z4 <- resp <- foo[, 2]
  z5 <- resp*(1-a1)*a2
  z6 <- resp*a1*(1-a2)
  y <- foo[, 4]
  fit_q2 <- lm(y ~ z1+z2+z3+z4+z5+z6)
  return (fit_q2)
}

FitPseudoOutcome <- function(q2, foo) {
  y <- foo[, 4]
  resp <- foo[, 2]
  n <- length(y)
  y_hat <- rep(NA, length(y))
  coef <- q2$coef
  coef[which(is.na(coef))] <- 0
  for (i in 1:n) {
    y_hat[i] <- Q2max(a1[i], resp[i], coef[1], coef[2], coef[3], coef[4], coef[5], coef[6], coef[7])[2]
  }
  return (y_hat)
}

FitQ1 <- function(y_hat, foo) {
  a1 <- foo[, 1]
  fit_q1 <- lm(y_hat ~ a1)
  return (fit_q1)
}

FitPi1 <- function(q1, foo, base) {
  a1 <- foo[, 1]
  q1_hat <- c(mean(q1$fit[a1 == 0]), mean(q1$fit[a1 == 1]))
  d1 <- q1_hat - min(q1_hat)
  d1 <- d1 / summary(q1)$sigma
  d1 <- pmin(d1, log(100000)/log(base))
  pi1 <- base^d1 / sum(base^d1)
  return (pi1)
}

FitPi2 <- function(q2, base) {
  Pi2 <- matrix(rep(NA, 16), nrow=4)
  colnames(Pi2) <- c("a1", "resp", "prob(a2=0)", "prob(a2=1)")
  coef <- q2$coef
  coef[which(is.na(coef))] <- 0
  iter <- 1
  for (a1 in 0:1) {
    for (r in 0:1) {
      q2hat <- c(Q2(a1, r, 0, coef[1], coef[2], coef[3], coef[4], coef[5], coef[6], coef[7]), Q2(a1, r, 1, coef[1], coef[2], coef[3], coef[4], coef[5], coef[6], coef[7]))
      d2 <- q2hat - min(q2hat)
      d2 <- d2 / summary(q2)$sigma
      d2 <- pmin(d2, log(100000)/log(base))
      prand2 <- base^d2 / sum(base^d2)
      Pi2[iter,] <- c(a1, r, prand2)
      iter <- iter + 1
    }
  }
  return (Pi2)
}


# Response Adaptive Randomisation ---------------------------
GetRandProb <- function(foo, base) {
  fit_q2 <- foo %>% FitQ2()
  fit_q1 <- fit_q2 %>%
    FitPseudoOutcome(foo = foo) %>%
    FitQ1(foo = foo)
  Pi1 <- fit_q1 %>% FitPi1(foo = foo, base = base)
  Pi2 <- fit_q2 %>% FitPi2(base = base)
  return (list(PI1 = Pi1, PI2 = Pi2))
}

GetTilProb <- function(w0, w1, pi0_at_0, pi0_at_1, pi_hat_at_0, pi_hat_at_1) {
  rho_0 <- exp(w0*log(pi0_at_0) + w1*log(pi_hat_at_0))
  rho_1 <- exp(w0*log(pi0_at_1) + w1*log(pi_hat_at_1))
  pi_til <- rho_1  / (rho_0 + rho_1)
  return (pi_til)
}


# Learning ---------------------------
Learning <- function(foo) {
  fit_q2 <- foo %>% FitQ2()
  y_hat <- fit_q2 %>% FitPseudoOutcome(foo = foo)
  fit_q1 <- y_hat %>% FitQ1(foo = foo)
  # values correspond to a1=0, a1=1
  q1_hat <- c(mean(fit_q1$fit[a1 == 0]), mean(fit_q1$fit[a1 == 1]))
  # values correspond to a1=0; r1=0,1
  #q20hat = c(Q2max(0,0, coef[1],coef[2],coef[3], coef[4],coef[5],coef[6],coef[7])[2], Q2max(0,1,coef[1],coef[2],coef[3],coef[4],coef[5],coef[6],coef[7])[2])
  # values correspond to a1=0; r1=0,1
  #q21hat = c(Q2max(1,0, coef[1],coef[2],coef[3], coef[4],coef[5],coef[6],coef[7])[2], Q2max(1,1,coef[1],coef[2],coef[3],coef[4],coef[5],coef[6],coef[7])[2])
  coef <- fit_q2$coef
  coef[which(is.na(coef))] <- 0

  odtr <- rep(NA, 3)
  if (q1_hat[1] >= q1_hat[2])  {
    odtr[1] <- 0
    odtr[2:3] <- c(Q2max(0,0, coef[1],coef[2],coef[3], coef[4],coef[5],coef[6],coef[7])[1], Q2max(0,1,coef[1],coef[2],coef[3],coef[4],coef[5],coef[6],coef[7])[1])
  }
  else {
    odtr[1] <- 1
    odtr[2:3] <- c(Q2max(1,0, coef[1],coef[2],coef[3], coef[4],coef[5],coef[6],coef[7])[1], Q2max(1,1,coef[1],coef[2],coef[3],coef[4],coef[5],coef[6],coef[7])[1])
  }
  return (list(optimal_dtr = odtr, fit_q1 = summary(fit_q1), fit_q2 = summary(fit_q2),q1_hat = q1_hat))
}


# Miscellaneous ---------------------------
GetLambda <- function(tau, b, N_min, N_complete) {
  lambda <- tau^{1/(b-1)} * N_min / N_complete
  return(lambda)
}
