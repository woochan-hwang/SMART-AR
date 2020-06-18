# Reimplementation and extension of SMART with adaptive randomisation [Biometrics, Cheung 2015]
# Woochan H. 2020, supervised by Dr David Robertson & Dr Sofia VIllar

# Load dependencies ------------------------
source("fit_models.R")
library(magrittr)

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
