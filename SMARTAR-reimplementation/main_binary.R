# Reimplementation and extension of SMART with adaptive randomisation [Biometrics, Cheung 2015]
# Woochan H. 2020, supervised by Dr David Robertson & Dr Sofia VIllar
# Following Google's R style guide

# LOAD DEPENDENCY ---------------------------
source("fit_models.R")
source("utility_functions.R")
library(magrittr)

# HYPERPARAMETERS ---------------------------
# Initial pi_0 fit to CODICAS-QoL data by Cheung
pi0_a1_1 <- pi0_a1_0 <- 0.5
pi0_a2_1 <- pi0_a2_0 <- c(0.5, 0.5, 0.5, 0.5)

# Patient arrival parameters
accrate <- 4   # patients / month
obswin <- 6  # trial duration (month), response analysed at (arrival + obswin)

# RAR algorithm design parameters
N_MIN <- 20
B <- 10
TAU <- 0.75

# Simulation parameters
N_PATIENTS <- 75
N_SIM <- 1000  # number of simulation repeats

FIXED_EQUAL_RANDOMISATION <- FALSE
if (FIXED_EQUAL_RANDOMISATION) {
  B <- 1
  pi0_a1_1 <- 0.5
  pi0_a1_0 <- 1 - pi0_a1_1
  pi0_a2_1 <- c(0.5, 0.5, 0.5, 0.5)
  pi0_a2_0 <- 1 - pi0_a2_1
}

# Probability of response to given action
response_prob_a1_a <- 0.5
response_prob_a1_b <- 0.5
response_prob_a2_aa <- 0.5
response_prob_a2_ab <- 0.4
response_prob_a2_ba <- 0.5
response_prob_a2_bb <- 0.5

# Hypothesis testing parameters
null_hypothesis <- "a2.given.a1"  # c("a1", "a2.combined", "a2.given.a1", "a2.given.r")
alternative <- "two.sided"  # c("two.sided", "greater", "less")
a1_h <- 0; r_h <- 0  # required for "a2.given.a1", "a2.given.r"
significance_level <- 0.1  # set as 0.1 in original paper

cat("[Session] True response rates: [a1]", response_prob_a1_a, response_prob_a1_b,
    "[a2]", response_prob_a2_aa, response_prob_a2_ab,
    response_prob_a2_ba, response_prob_a2_bb, "[null]", null_hypothesis, "\n")

# True outcome based on given parameters
ground_truth <- ExperimentSetting(null_hypothesis, response_prob_a1_a, response_prob_a1_b, response_prob_a2_aa,
                                  response_prob_a2_ab, response_prob_a2_ba, response_prob_a2_bb, a1_h)$ground_truth
superior_treatment <- ExperimentSetting(null_hypothesis, response_prob_a1_a, response_prob_a1_b, response_prob_a2_aa,
                                        response_prob_a2_ab, response_prob_a2_ba, response_prob_a2_bb, a1_h)$superior_treatment


# RUN SIMULATIONS ------------------------------
set.seed(0202)
pval <- patient_outcome <- c(rep(NA, N_SIM))
cat("[Running]", N_SIM, "simulation iterations for", N_PATIENTS, "samples \n")
ptm <- proc.time()

for (sim in 1:N_SIM) {

  arrival <- c(0, cumsum(rexp(N_PATIENTS - 1, rate = accrate)))
  a1 <- a2 <- R <- y <- rep(NA, N_PATIENTS)

  for (i in 1:N_PATIENTS) {
    ind_comp <- which(arrival < (arrival[i] - obswin))  # True when patient outcome is available (arrival + trial duration)
    n_comp <- length(ind_comp)

    if (n_comp < N_MIN)  {
      # Generate DTR sequence
      a1[i] <- Action1(pi1 = pi0_a1_1)
      R[i] <- InterimResponse(a1[i], p0 = response_prob_a1_a, p1 = response_prob_a1_b)
      a2[i] <- Action2(a1[i], R[i], pi2 = pi0_a2_1)
      y[i] <- GetOutcome(a1[i], a2[i], response_prob_a2_aa, response_prob_a2_ab, response_prob_a2_ba, response_prob_a2_bb)
    }
    else {
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
      R[i] <- InterimResponse(a1[i], p0 = response_prob_a1_a, p1 = response_prob_a1_b)
      a2[i] <- Action2(a1[i], R[i], pi2 = pi2_til)
      y[i] <- GetOutcome(a1[i], a2[i], response_prob_a2_aa, response_prob_a2_ab, response_prob_a2_ba, response_prob_a2_bb)
    }
  }

  df <- data.frame(a1, R, a2, y)
  #PrintIterationSummary(df)
  pval[sim] <- HypothesisTesting(dataframe = df, null = null_hypothesis, alternative = alternative,
                                 a1.h = a1_h, r.h = r_h)$p.value
  patient_outcome[sim] <- PatientOutcome(dataframe = df, null = null_hypothesis, a1.h = a1_h, r.h = r_h,
                                         superior.treatment = superior_treatment)
}

# PRINT OUTCOME ---------------------------
cat("[Completed]", N_SIM, "simulation iterations:", (proc.time() - ptm)[[3]], "s \n")

# HYPOTHESIS TESTING ------------------------
if (ground_truth) {
  type_1_error <- length(which(pval < significance_level)) / N_SIM
  cat("Type 1 Error:", type_1_error, "\n")
} else {
  power <- length(which(pval < significance_level)) / N_SIM
  cat("Statistical Power:", power, "\n")
}

cat("Proportion of patients allocated to superior arm:", mean(patient_outcome), "\n")

