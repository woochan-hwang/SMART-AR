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
N_PATIENTS <- 100
N_SIM <- 2  # number of simulation repeats


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
    print(ind_comp)

    if (n_comp < N_MIN)  {
      # Generate DTR sequence
      a1[i] <- Action1(pi1 = pi0_a1_1)
      R[i] <- InterimResponse(a1[i], p0 = response_prob_a1_0, p1 = response_prob_a1_1)
      a2[i] <- Action2(a1[i], R[i], pi2 = pi0_a2_1)
    }
    else {
      y <- Q2sat(a1,R,a2, theta_sat) + epsilon
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

      pi_a1 <- c(pi0_a1_0, pi0_a1_1, pi_hat_a1_0, pi_hat_a1_1)
      pi1_til <- GetTilProb(w0 = w0, w1 = w1, pi_a1)

      pi_a2 <- c(pi0_a2_0, pi0_a2_1, pi_hat_a2_0, pi_hat_a2_1)
      pi2_til <- GetTilProb(w0 = w0, w1 = w1, pi_a2)

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
cat("\nSMART-AR continuous monitoring\n")
cat("N:\t",N_PATIENTS,"\n")
cat("Nmin:",N_MIN,"\n")
cat("base:",B,"\n")
cat("tau:",TAU,"\n")
cat("pi1:",pi0_a1_1,"\n")
cat("pi2:",pi0_a2_1,"\n")
cat("arate:", accrate,"\n")
cat("nsim:",N_SIM,"\n\n")

v4 <- dtr_table[, 4]
vmax <- which(v4 == max(v4))
vmax <- dtr_table[vmax, 4]
vmin <- which(v4 == min(v4))
vmin <- dtr_table[vmin, 4]
pcs <- length(which(abs(val-vmax) <= 0.0001))/N_SIM
pcsg <- length(which(abs(valg-vmax) <= 0.0001))/N_SIM

cat("Distribution of the values of the selected DTR (dhat):\n")
value_dhat <- val
print(summary(value_dhat))
print(table(value_dhat),digits=3)
cat("\n")

cat("Distribution of the average patient outcome per trial\n")
patientval <- apply(Y, 1, mean)
print(summary(patientval))

cat("Some key summary statistics\n")
cat("\nProbability of selecting the true optimal DTR:", pcs, "\n\n")
cat("Value of selected dtr (Mean):", mean(val),"\n")
cat("Value of selected dtr (SD):", sd(val),"\n")
cat("Value of selected dtr (Q1,Q2,Q3)", quantile(val,c(0.25,0.5,0.75)),"\n")
cat("Adjusted value of selected dtr (Mean)", (mean(val)-vmin)/(vmax-vmin),"\n")
cat("Adjusted value of selected dtr (SD)", sd(val)/(vmax-vmin),"\n")
cat("Adjusted value of selected dtr (Q1,Q2,Q3)", quantile(val-vmin,c(0.25,0.5,0.75))/(vmax-vmin),"\n\n")

cat("Average patient outcome (Mean):", mean(patientval),"\n")
cat("Average patient outcome (SD):", sd(patientval),"\n")
cat("Average patient outcome (Q1,Q2,Q3):", quantile(patientval,c(0.25,0.5,0.75)),"\n")
cat("Adjusted average patient outcome (Mean):", (mean(patientval)-vmin)/(vmax-vmin),"\n")
cat("Adjusted average patient outcome (SD):", sd(patientval)/(vmax-vmin),"\n")
cat("Adjusted average patient outcome (Q1,Q2,Q3):", quantile(patientval-vmin,c(0.25,0.5,0.75))/(vmax-vmin),"\n")
