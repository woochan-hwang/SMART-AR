# Reimplementation and extension of SMART with adaptive randomisation [Biometrics, Cheung 2015]
# Woochan H. 2020, supervised by Dr David Robertson & Dr Sofia VIllar
# Following Google's R style guide

# LOAD DEPENDENCY ---------------------------
source("fit_models.R")
source("utility_functions.R")
library(ggplot2)
library(magrittr)

# HYPERPARAMETERS ---------------------------
# Initial pi_0 fit to CODICAS-QoL data by Cheung
pi0_a1_1 <- 0.6749323
pi0_a1_0 <- 1 - pi0_a1_1
pi0_a2_1 <- c(0.7016099, 0.3799304, 0.4030912, 0.2576907)
pi0_a2_0 <- 1 - pi0_a2_1

# Probability of response to a1 in CODIACS-QoL
response_prob_a1_0 <- 29 / 56
response_prob_a1_1 <- 28 / 52
outcome_sigma <- 6.7

# Patient arrival parameters
accrate <- 4   # patients / month
obswin <- 6  # trial duration (month), response analysed at (arrival + obswin)

# Regression parameters used by Cheung
THETA_1 <- c(2.206, 5.594, 8.294, 7.745, -12.102, 6.455, -13.045, 6.59)     # Scene 1
THETA_2 <- c(2.206, 5.594, 8.294, 7.745, -6.051, 6.455, -6.5225, 0.0675)    # Scene 2
THETA_3 <- c(2.206, 5.594, 8.294, 7.745, -12.102, -6.455, -13.045, 19.5)    # Scene 3
THETA_4 <- c(2.206, 5.594, 8.294, 7.745, -6.051, -6.455, -6.5225, 12.9775)  # Scene 4
THETA_5 <- c(1.320, 6.480, 9.180, 9.555, -11.822, 4.645, -14.855, 6.382)    # Scene 5 (fitted saturated model)
THETA_6 <- c(2.206, 5.594, 8.294, 7.745, -12.102, -6.455, -13.045, -19.5)    # Scene 6 (modified scene 3)

# RAR algorithm design parameters
N_MIN <- 20
B <- 10
TAU <- 0.75

# Simulation parameters
N_PATIENTS <- 75
N_SIM <- 1  # number of simulation repeats

FIXED_EQUAL_RANDOMISATION <- FALSE
if (FIXED_EQUAL_RANDOMISATION) {
  B <- 1
  pi0_a1_1 <- 0.5
  pi0_a1_0 <- 1 - pi0_a1_1
  pi0_a2_1 <- c(0.5, 0.5, 0.5, 0.5)
  pi0_a2_0 <- 1 - pi0_a2_1
}

# PRINT STARTING POINT -------------------------
theta_sat <- THETA_3
dtr_table <- DtrTable(response_prob_a1_0, response_prob_a1_1, theta_sat) %>% PrintDtrTable()

# RUN SIMULATIONS ------------------------------
set.seed(0202)
DTRg <- DTRhat <- matrix(rep(NA, N_SIM * 3), nrow = N_SIM)
valg <- val <- rep(NA, N_SIM)
Y <- matrix(rep(NA, N_SIM * N_PATIENTS), nrow = N_SIM)

for (r in 1:N_SIM) {
  if (round(r / 100) == (r / 100)) print(r)

  # Patient outcomes and arrival times
  epsilon <- rnorm(N_PATIENTS, 0, outcome_sigma)
  arrival <- c(0, cumsum(rexp(N_PATIENTS - 1, rate = accrate)))

  # Sequential experiment (a1, a2)
  a1 <- a2 <- R <- y <- rep(NA, N_PATIENTS)
  for (i in 1:N_PATIENTS) {
    ind_comp <- which(arrival < (arrival[i] - obswin))  # True when patient outcome is available (arrival + trial duration)
    n_comp <- length(ind_comp)

    if (n_comp < N_MIN)  {
      # Generate DTR sequence
      a1[i] <- Action1(pi1 = pi0_a1_1)
      R[i] <- InterimResponse(a1[i], p0 = response_prob_a1_0, p1 = response_prob_a1_1)
      a2[i] <- Action2(a1[i], R[i], pi2 = pi0_a2_1)
    }
    else {
      y <- Q2sat(a1,R,a2, theta_sat) + epsilon
      cat("y", y)
      a1_comp <- a1[ind_comp]
      R_comp <- R[ind_comp]
      a2_comp <- a2[ind_comp]
      y_comp <- y[ind_comp]

      # regression model to calculate pi
      foo <- cbind(a1_comp,R_comp,a2_comp,y_comp)
      pfoo <- GetRandProb(foo, base = B)
      pi_hat_a1_0 <- pfoo$PI1[1]
      pi_hat_a1_1 <- pfoo$PI1[2]
      pi_hat_a2_0 <- pfoo$PI2[, 3]
      pi_hat_a2_1 <- pfoo$PI2[, 4]

      # weighted average for historical context integration
      lambda <- GetLambda(tau = TAU, b = B, N_min = N_MIN, N_complete = n_comp)
      w0 <- lambda^{B-1}
      w1 <- 1 - w0

      pi1_til <- GetTilProb(w0 = w0, w1 = w1, pi0_a1_0, pi0_a1_1, pi_hat_a1_0, pi_hat_a1_1)
      pi2_til <- GetTilProb(w0 = w0, w1 = w1, pi0_a2_0, pi0_a2_1, pi_hat_a2_0, pi_hat_a2_1)

      # Generate DTR sequence
      a1[i] <- Action1(pi1 = pi1_til)
      R[i] <- InterimResponse(a1[i], p0 = response_prob_a1_0, p1 = response_prob_a1_1)
      a2[i] <- Action2(a1[i], R[i], pi2 = pi2_til)
    }
  }

  # Final fit after complete follow-up of all subjects
  y <- Q2sat(a1, R, a2, theta_sat) + epsilon
  foo <- cbind(a1, R, a2, y)
  optimal_fit <- Learning(foo)

  Y[r,] <- y
  DTRhat[r,] <- dtrhat <- optimal_fit$optimal_dtr

  pos <- which(dtr_table[, 1] == dtrhat[1] &
                 dtr_table[, 2] == dtrhat[2] &
                 dtr_table[, 3] == dtrhat[3])
  val[r] <- dtr_table[pos, 4]  #optimal dtr Q-value
  #DTRg[r,] <- GetGEstimate(foo)$odtr
  #pos <- which(dtr_table[, 1] == DTRg[r, 1] &
  #               dtr_table[, 2] == DTRg[r, 2] &
  #               dtr_table[, 3] == DTRg[r, 3])
  #valg[r] = dtr_table[pos, 4]
}

# PRINT OUTCOME ---------------------------
PrintHyperParameters(theta_sat, N_MIN, B, TAU, N_PATIENTS, N_SIM, pi0_a1_1, pi0_a2_1, accrate, obswin)
PrintSummary(dtr_table, val)