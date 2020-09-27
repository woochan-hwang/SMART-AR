# Applying constrained randomised dynamic programming to SMART
# Woochan H. 2020, supervised by Dr David Robertson & Dr Sofia VIllar

# LOAD DEPENDENCY ---------------------------
suppressMessages(library(R.utils))
source("nested_crdp.R")
source("utility_functions.R")
options(warnings=-1)


# HYPERPARAMETERS ---------------------------
# RAR algorithm design parameters
degree_of_randomisation <- 0.9
degree_of_constraint <- 0.15  # DoC * N_patients / N_trial_arms = minimum allocation per DTR

# Prior probabilities
bayes_prior_a1 <- c(1,1,1,1)  # prior used for crdp value function
bayes_prior_a2 <- c(1,1,1,1)

# Simulation parameters
N_SIM <- 100  # number of simulation repeats
N_PATIENTS <- 10

# True distribution parameters
response_prob_a1_a <- 0.7
response_prob_a1_b <- 0.5
response_prob_a2_aa <- 0.5
response_prob_a2_ab <- 0.1
response_prob_a2_ba <- 0.5
response_prob_a2_bb <- 0.5

# Hypothesis testing parameters
null_hypothesis <- "a2.given.a1"  # c("a1", "a2.combined", "a2.given.a1", "a2.given.r")
alternative <- "two.sided"  # c("two.sided", "greater", "less")
a1_h <- 0; r_h <- 0  # required for "a2.given.a1", "a2.given.r"
significance_level <- 0.1  # set as 0.1 in original paper

# Overwrite based on command line input
args <- commandArgs(asValue=TRUE, excludeReserved=TRUE, adhoc=TRUE)[-1]
keys <- attachLocally(args, overwrite=TRUE)

# True outcome based on given parameters
ground_truth <- ExperimentSetting(null_hypothesis, response_prob_a1_a, response_prob_a1_b, response_prob_a2_aa,
                                  response_prob_a2_ab, response_prob_a2_ba, response_prob_a2_bb, a1_h)$ground_truth
superior_treatment <- ExperimentSetting(null_hypothesis, response_prob_a1_a, response_prob_a1_b, response_prob_a2_aa,
                                        response_prob_a2_ab, response_prob_a2_ba, response_prob_a2_bb, a1_h)$superior_treatment
parameterList <- PrintHyperParameters(returnList=TRUE)


# RUN SIMULATIONS ------------------------------
set.seed(0202)
pval <- patient_outcome <- c(rep(NA, N_SIM))
# calculate optimal policy via nested CRDP
# input: (sample_size, a1_prior_vector, a2_prior_vector, randomisation, constraint) where the prior vectors
# are vectors of length 4 (s_a,f_a,s_b,f_b); 0 <= constraint <= 0.5 <= randomisation <= 1; sample_size = int;
# return: a list with "Q1" and "Q2", where Q1 <- list("ACTION", "VALUE"); Q2 <- list("N_i_ACTION") for i in range(Q2)
optimal_policy <- TrainNestedCRDP(sample_size = N_PATIENTS, a1_prior_vector = bayes_prior_a1, a2_prior_vector = bayes_prior_a2,
                                  randomisation = degree_of_randomisation, constraint = degree_of_constraint)

cat("[Running]", N_SIM, "simulation iterations for", N_PATIENTS, "samples \n")
ptm <- proc.time()

for (sim in 1:N_SIM) {

  i <- j <- k <- l <- 1
  a1 <- r <- a2 <- y <- rep(NA, N_PATIENTS)

  # DTR Step 1 Inference
  for (t in 1:N_PATIENTS) {
    n <- N_PATIENTS - t + 1
    optimal_action <- optimal_policy$Q1$ACTION[n, i, j, k]  # optimal_action <- 0 if "treatment_a"; 1 if "treatment_b"
    a1[t] <- selected_action <- SelectTreatment(optimal_action, degree_of_randomisation)
    r[t] <- GetInterimResponse(selected_action, response_prob_a1_a, response_prob_a1_b)
    belief_a1 <- UpdateCurrentBelief(c(i,j,k,l), a1[t], r[t])
    i <- belief_a1[1]; j <- belief_a1[2]; k <- belief_a1[3]; l <- belief_a1[4]
  }

  O2_group_size <- belief_a1 - 1
  patient_index_list <- list(which(a1==0)[which(r[a1==0]==1)], which(a1==0)[which(r[a1==0]==0)],
                             which(a1==1)[which(r[a1==1]==1)], which(a1==1)[which(r[a1==1]==0)])

  # DTR Step 2 Inference
  for (group in 1:4) {
    group_size <- O2_group_size[group]
    if (group_size == 0) {next}  # skip iteration for groups with no patients
    group_name <- paste("N", group_size, "ACTION", sep="_")
    if (group_size < 4) {group_name <- "N_4_ACTION"}
    optimal_group_policy <- optimal_policy$Q2[[group_name]]
    i <- j <- k <- l <- 1

    # A2 selection based on optimal policy per subgropu
    for (t in 1:group_size) {
      patient_index <- patient_index_list[[group]][t]
      n <- group_size - t + 1
      optimal_action <- optimal_group_policy[n, i, j, k]
      a2[patient_index] <- selected_action <- SelectTreatment(optimal_action, degree_of_randomisation)
      y[patient_index] <- GetOutcome(a1[patient_index], selected_action, response_prob_a2_aa, response_prob_a2_ab,
                                     response_prob_a2_ba, response_prob_a2_bb)
      belief_a2 <- UpdateCurrentBelief(c(i,j,k,l), a2[patient_index], y[patient_index])
      i <- belief_a2[1]; j <- belief_a2[2]; k <- belief_a2[3]; l <- belief_a2[4]
    }
  }

  df <- data.frame(a1, r, a2, y)
  #PrintIterationSummary(df)
  pval[sim] <- HypothesisTesting(dataframe = df, null = null_hypothesis, alternative = alternative,
                                 a1.h = a1_h, r.h = r_h)$p.value
  patient_outcome[sim] <- PatientOutcome(dataframe = df, null = null_hypothesis, a1.h = a1_h, r.h = r_h,
                                         superior.treatment = superior_treatment)
}

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
