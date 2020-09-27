# Reimplementation and extension of SMART with adaptive randomisation [Biometrics, Cheung 2015]
# Woochan H. 2020, supervised by Dr David Robertson & Dr Sofia VIllar

# Load dependencies ------------------------
source("fit_models.R")
library(magrittr)
library(assertthat)

# DTR logic functions ----------------------
Action1 <- function(pi1) {
  return (rbinom(1, 1, pi1))
}

InterimResponse <- function(a1, p0, p1) {
  if (a1 == 0) {
    R <- rbinom(1, 1, p0)
  } else {
    R <- rbinom(1,1,p1)
  }
  return (R)
}

Action2 <- function(a1, R, pi2) {
  # STAGE 2 of DTR
  if (a1 == 0 & R == 0) {
    a2 <- rbinom(1,1,pi2[1])
  } else if (a1 == 0 & R == 1) {
    a2 <- rbinom(1,1,pi2[2])
  } else if (a1 == 1 & R == 0) {
    a2 <- rbinom(1,1,pi2[3])
  } else {
    a2 <- rbinom(1,1,pi2[4])
  }
  return (a2)
}


# Hypothesis testing -------------------------

# standard CI construction for a sample drawn from a normal distribution
ConfidenceInterval <- function(sample, range = 95) {
  mean <- mean(sample)
  standard_dev <- sd(sample)
  sample_size <- length(sample)
  z <- c(2.576, 2.326, 1.96, 1.645)[which(range == c(99, 98, 95, 90))]
  lower_bound <- mean - z * standard_dev / sqrt(sample_size)
  upper_bound <- mean + z * standard_dev / sqrt(sample_size)
  return (c(lower_bound, upper_bound))
}


# Output functions ---------------------------
DtrTable <- function(response_prob_a1_0, response_prob_a1_1, theta_sat) {
  dtr_table <- matrix(rep(NA,4*8),nrow=8)
  dtr_table[,1] <- c(rep(0,4),rep(1,4))
  dtr_table[,2] <- c(0,0,1,1,0,0,1,1)
  dtr_table[,3] <- c(0,1,0,1,0,1,0,1)

  for (i in 1:8) {
    a1 <- dtr_table[i,1]
    if (a1 == 0) {
      r1 <- response_prob_a1_0
    } else if (a1 == 1) {
      r1 <- response_prob_a1_1
    }
    r0 <- 1 - r1
    dtr_table[i,4] <- r0*Q2sat(a1,0,dtr_table[i,2], theta_sat) + r1*Q2sat(a1,1,dtr_table[i,3], theta_sat)
  }
  return (dtr_table)
}

PrintDtrTable <- function(dtr_table) {
  dtr_value <- dtr_table[,4]
  best_pos <- which(dtr_value==max(dtr_value))
  worst_pos <- which(dtr_value==min(dtr_value))
  colnames(dtr_table) <- c("d1","d2(0)","d2(1)","Value")
  row <- rep("", nrow(dtr_table))
  row[best_pos] <- "OPTIMAL"
  row[worst_pos] <- "WORST"
  rownames(dtr_table) <- row

  cat("\nAll possible DTRs and their values","\n")
  print(dtr_table, digits=3)
  invisible(dtr_table)
}

PrintHyperParameters <- function(theta_sat, N_MIN, B, TAU, N_PATIENTS, N_SIM, pi0_a1_1, pi0_a2_1, accrate, obswin) {
  # print regression params
  cat("REGRESSION PARAMETERS OF Saturated model:\n")
  thetaprn <- matrix(theta_sat, nrow=1)
  colnames(thetaprn) <- c("beta0", "beta1", "beta2", "gamma1", "beta3", "gamma3", "gamma2", "gamma4")
  print(thetaprn)
  # print other hyperparameters
  cat("\nSMART-AR patient parameters\n")
  cat("N:\t", N_PATIENTS, "\n")
  cat("arate:\t", accrate, "\n")
  cat("obswin:\t", obswin, "\n")
  cat("pi1:\t", pi0_a1_1, "\n")
  cat("pi2:\t", pi0_a2_1, "\n")

  cat("\nSMART-AR simulation parameters\n")
  cat("N_sim:\t", N_SIM, "\n")
  cat("Nmin:\t", N_MIN, "\n")
  cat("base:\t", B, "\n")
  cat("tau:\t", TAU, "\n\n")
}

PrintSummary <- function(dtr_table, val) {
  v4 <- dtr_table[, 4]
  vmax <- which(v4 == max(v4))
  vmax <- dtr_table[vmax, 4]
  vmin <- which(v4 == min(v4))
  vmin <- dtr_table[vmin, 4]
  pcs <- length(which(abs(val-vmax) <= 0.0001))/N_SIM

  cat("Distribution of the values of the selected DTR (dhat):\n")
  print(summary(val))
  print(table(val),digits=3)
  cat("\n")

  cat("Distribution of the average patient outcome per trial\n")
  patientval <- apply(Y, 1, mean)
  print(summary(patientval))

  cat("Some key summary statistics\n")
  cat("\nProbability of selecting the true optimal DTR:", pcs, "\n\n")
  cat("Value of selected dtr (Mean):", mean(val),"\n")
  cat("Value of selected dtr (SD):", sd(val),"\n")
  cat("Value of selected dtr (Q1,Q2,Q3)", quantile(val,c(0.25,0.5,0.75)),"\n")
  cat("95% Confidence Interval for selected dtr:", ConfidenceInterval(val, 95),"\n")
  cat("Adjusted value of selected dtr (Mean)", (mean(val)-vmin)/(vmax-vmin),"\n")
  cat("Adjusted value of selected dtr (SD)", sd(val)/(vmax-vmin),"\n")
  cat("Adjusted value of selected dtr (Q1,Q2,Q3)", quantile(val-vmin,c(0.25,0.5,0.75))/(vmax-vmin),"\n\n")

  cat("Average patient outcome (Mean):", mean(patientval),"\n")
  cat("Average patient outcome (SD):", sd(patientval),"\n")
  cat("Average patient outcome (Q1,Q2,Q3):", quantile(patientval,c(0.25,0.5,0.75)),"\n")
  cat("Adjusted average patient outcome (Mean):", (mean(patientval)-vmin)/(vmax-vmin),"\n")
  cat("Adjusted average patient outcome (SD):", sd(patientval)/(vmax-vmin),"\n")
  cat("Adjusted average patient outcome (Q1,Q2,Q3):", quantile(patientval-vmin,c(0.25,0.5,0.75))/(vmax-vmin),"\n")
}

# PERFORMANCE METRICS -----------------------
# Fisher's Exact Test, small samples of two binomial distributions
HypothesisTesting <- function(dataframe, null, alternative, a1.h = NA, r.h = NA) {

  if (null == "a2.given.a1") {assert_that(!is.na(a1.h))}
  if (null == "a2.given.r") {assert_that(!is.na(a1.h), !is.na(r.h))}
  a1 <- dataframe$a1
  r <- dataframe$R
  a2 <- dataframe$a2
  y <- dataframe$y

  if (null == "a1") {
    sim_results <- rbind(c(length(which(r[a1==0]==1)), length(which(r[a1==0]==0))),
                         c(length(which(r[a1==1]==1)), length(which(r[a1==1]==0))))
  } else if (null == 'a2.combined') {
    sim_results <- rbind(c(length(which(y[a2==0]==1)), length(which(y[a2==0]==0))),
                         c(length(which(y[a2==1]==1)), length(which(y[a2==1]==0))))
  } else if (null == 'a2.given.a1') {
    sim_results <- rbind(c(length(which(y[which(a2[a1==a1.h]==0)]==1)), length(which(y[which(a2[a1==a1.h]==0)]==0))),
                         c(length(which(y[which(a2[a1==a1.h]==1)]==1)), length(which(y[which(a2[a1==a1.h]==1)]==0))))
  } else if (null == 'a2.given.r') {
  }
  #  cat("Hypothesis testing performed on: ", sim_results, "\n")
  results <- fisher.test(sim_results, alternative = alternative, conf.int = FALSE)
  return (results)
}

PatientOutcome <- function(dataframe, null, a1.h, r.h, superior.treatment) {

  if (null == "a2.given.a1") {assert_that(!is.na(a1.h))}
  a1 <- dataframe$a1
  a2 <- dataframe$a2

  if (null == "a1") {
    outcome <- length(which(a1==superior.treatment)) / length(a1)
  } else if (null == 'a2.combined') {
    outcome <- length(which(a2==superior.treatment)) / length(a2)
  } else if (null == 'a2.given.a1') {
    outcome <- length(which(a2[a1==a1.h]==superior.treatment)) / length(which(a1==a1.h))
  } else if (null == 'a2.given.r') {
  }
  return (outcome)
}

ExperimentSetting <- function(null, response_prob_a1_a, response_prob_a1_b, response_prob_a2_aa, response_prob_a2_ab, response_prob_a2_ba, response_prob_a2_bb, a1_h) {

  if (a1_h == 0) {a1_h <- 'a'}
  else {a1_h <- 'b'}

  if (null == "a1") {
    prob1 <- response_prob_a1_a
    prob2 <- response_prob_a1_b
  } else if (null == 'a2.combined') {
    prob1 <- (response_prob_a2_aa + response_prob_a2_ba) / 2
    prob2 <- (response_prob_a2_ab + response_prob_a2_bb) / 2
  } else if (null == 'a2.given.a1') {
    prob1 <- get(paste("response_prob_a2_", a1_h, "a", sep=""))
    prob2 <- get(paste("response_prob_a2_", a1_h, "b", sep=""))
  } else if (null == 'a2.given.r') {
    # look up interaction design
  }

  if (prob1 == prob2) {ground_truth <- TRUE} else {ground_truth <- FALSE}
  if (prob1 > prob2) {superior_treatment <- 0} else {superior_treatment <- 1}

  return(list("ground_truth"=ground_truth, "superior_treatment"=superior_treatment))
}


# KEY INFERENCE FUNCTIONS ---------------------------
SelectTreatment <- function(optimal_action, degree_of_randomisation) {
  if (optimal_action == 0) {
    selected_action <- 1 - rbinom(1, 1, degree_of_randomisation)  # d_o_r = P(opitmal_action == selected_action)
  }
  else if (optimal_action == 1) {
    selected_action <- rbinom(1,1, degree_of_randomisation)
  }
  else {
    selected_action <- rbinom(1,1, 0.5)  # optimal_action == 2 when ambiguous
  }
  return (selected_action)
}

GetInterimResponse <- function(selected_action, response_prob_a1_a, response_prob_a1_b) {
  if (selected_action == 0) {
    response <- rbinom(1, 1, response_prob_a1_a)
  }
  if (selected_action == 1) {
    response <- rbinom(1, 1, response_prob_a1_b)
  }
  return (response)
}

UpdateCurrentBelief <- function(belief_vector, action, response) {
  i <- belief_vector[1]
  j <- belief_vector[2]
  k <- belief_vector[3]
  l <- belief_vector[4]

  if (action == 0) {
    if (response == 1) {i <- i + 1}
    if (response == 0) {j <- j + 1}
  }
  if (action == 1) {
    if (response == 1) {k <- k + 1}
    if (response == 0) {l <- l + 1}
  }
  return(c(i,j,k,l))
}

GetOutcome <- function(a1, a2, response_prob_a2_aa, response_prob_a2_ab, response_prob_a2_ba, response_prob_a2_bb) {
  if (a2 == 0) {
    # generate response given (o1, a1, o2, a2)
    if (a1 == 0) {response <- rbinom(1, 1, response_prob_a2_aa)}
    if (a1 == 1) {response <- rbinom(1, 1, response_prob_a2_ba)}
  }
  if (a2 == 1) {
    if (a1 == 0) {response <- rbinom(1, 1, response_prob_a2_ab)}
    if (a1 == 1) {response <- rbinom(1, 1, response_prob_a2_bb)}
  }
  return (response)
}