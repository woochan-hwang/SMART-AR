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
