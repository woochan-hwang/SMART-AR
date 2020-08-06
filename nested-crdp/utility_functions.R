# Applying constrained randomised dynamic programming to SMART
# Woochan H. 2020, supervised by Dr David Robertson & Dr Sofia VIllar
# Includes misc functions for hypothesis testing of simulation outputs

# LOAD DEPENDENCY ---------------------------


# PERFORMANCE METRICS -----------------------
# Fisher's Exact Test, small samples of two binomial distributions
HypothesisTesting <- function(dataframe, null, alternative, a1.h = NA, r.h = NA) {

  if (null == "a2.given.a1") {assert_that(!is.na(a1.h))}
  if (null == "a2.given.r") {assert_that(!is.na(a1.h), !is.na(r.h))}
  a1 <- dataframe$a1
  r <- dataframe$r
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
    prob1 <- paste("response_prob_a2_", a1_h, "a", sep="")
    prob2 <- paste("response_prob_a2_", a1_h, "b", sep="")
  } else if (null == 'a2.given.r') {
    # look up interaction design
  }

  if (prob1 == prob2) {ground_truth <- TRUE} else {ground_truth <- FALSE}
  if (prob1 > prob2) {superior_treatment <- 0} else {superior_treatment <- 1}

  return(list("ground_truth"=ground_truth, "superior_treatment"=superior_treatment))
}


# OUTPUT FUNCTIONS --------------------------
PrintIterationSummary <- function(dataframe) {
  names(dataframe) <- c("a1", "r", "a2", "y")
  print(head(dataframe, 5))
  a1 <- dataframe$a1
  r <- dataframe$r

  cat("DTR step 1 Summary: \n")
  a1_outcome <- list(which(a1==0)[which(r[a1==0]==1)], which(a1==0)[which(r[a1==0]==0)],
                     which(a1==1)[which(r[a1==1]==1)], which(a1==1)[which(r[a1==1]==0)])
  cat("A1 | treatment a =", length(which(a1==0)), "; treatment b =", length(which(a1==1)), "\n")
  cat("O2 | treatment a = (", length(a1_outcome[[1]]), ",", length(a1_outcome[[2]]),
      "); treatment b = (", length(a1_outcome[[3]]), ",", length(a1_outcome[[4]]), ") \n")
  cat("DTR step 2 Summary: \n")
  a2_outcome <- list("(a,1,a)"=which(a2==0)[which(a2[a1_outcome[[1]]]==0)], "(a,1,b)"=which(a2==1)[which(a2[a1_outcome[[1]]]==1)],
                     "(a,0,a)"=which(a2==0)[which(a2[a1_outcome[[2]]]==0)], "(a,0,b)"=which(a2==1)[which(a2[a1_outcome[[2]]]==1)],
                     "(b,1,a)"=which(a2==0)[which(a2[a1_outcome[[3]]]==0)], "(b,1,b)"=which(a2==1)[which(a2[a1_outcome[[3]]]==1)],
                     "(b,0,a)"=which(a2==0)[which(a2[a1_outcome[[4]]]==0)], "(b,0,b)"=which(a2==1)[which(a2[a1_outcome[[4]]]==1)])
#  cat("A2 |", "\n")
#  print(a2_outcome)
#  cat("\n")
#  cat("Y  |", y, "\n")
  return(dataframe)
}

PrintHyperParameters <- function(returnList) {
  listOfParameters <- list()
  parameterNames <- c("degree_of_randomisation", "degree_of_constraint", "bayes_prior_a1", "bayes_prior_a2",
                      "N_SIM", "N_PATIENTS", "response_prob_a1_a", "response_prob_a1_b", "response_prob_a2_aa",
                      "response_prob_a2_ab", "response_prob_a2_ba", "response_prob_a2_bb", "null_hypothesis",
                      "alternative", "a1_h", "r_h", "significance_level", "ground_truth", "superior_treatment")
  for (obj in parameterNames) {
    listOfParameters[[obj]] <- get(obj)
  }
  listOfParameters <- as.data.frame(listOfParameters)
  print(head(listOfParameters, 1))
  if (returnList) {return (listOfParameters)}
}

SaveExperimentOutcome <- funciton(parameterList, overwrite=FALSE) {

}