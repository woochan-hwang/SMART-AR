# Applying constrained randomised dynamic programming to SMART
# Woochan H. 2020, supervised by Dr David Robertson & Dr Sofia VIllar

# LOAD DEPENDENCY ---------------------------
source("nested_crdp.R")
library(magrittr)

# HYPERPARAMETERS ---------------------------
# RAR algorithm design parameters
degree_of_randomisation <- 0.7
degree_of_constraint <- 0.3  # DoC * N_patients / N_trial_arms = minimum allocation per DTR

# Prior probabilities
bayes_prior_a1 <- c(1,1,1,1)  # prior used for crdp value function
bayes_prior_a2 <- c(1,1,1,1)

# Simulation parameters
N_SIM <- 1  # number of simulation repeats
N_PATIENTS <- 20
#accrate <- 4   # patients / month
#obswin <- 6  # trial duration (month), response analysed at (arrival + obswin)

# True distribution parameters
response_prob_a1_a <- 0.5
response_prob_a1_b <- 0.5
response_prob_a2_aa <- 0.4
response_prob_a2_ab <- 0.8
response_prob_a2_ba <- 0.4
response_prob_a2_bb <- 0.4


# RUN SIMULATIONS ------------------------------
set.seed(0202)
# calculate optimal policy via nested CRDP
# input: (sample_size, a1_prior_vector, a2_prior_vector, randomisation, constraint) where the prior vectors
# are vectors of length 4 (s_a,f_a,s_b,f_b); 0 <= constraint <= 0.5 <= randomisation <= 1; sample_size = int;
# return: a list with "Q1" and "Q2", where Q1 <- list("ACTION", "VALUE"); Q2 <- list("N_i_ACTION") for i in range(Q2)
optimal_policy <- TrainNestedCRDP(sample_size = N_PATIENTS, a1_prior_vector = bayes_prior_a1, a2_prior_vector = bayes_prior_a2,
                                  randomisation = degree_of_randomisation, constraint = degree_of_constraint)

for (sim in 1:N_SIM) {

  a1 <- r <- a2 <- y <- rep(NA, N_PATIENTS)
  i <- j <- k <- l <- 1

  # DTR Step 1 Inference
  for (t in 1:N_PATIENTS) {
    n <- N_PATIENTS - t + 1
    optimal_action <- optimal_policy$Q1$ACTION[n, i, j, k]  # optimal_action <- 0 if "treatment_a"; 1 if "treatment_b"

    a1[t] <- selected_action <- SelectTreatment(optimal_action, degree_of_randomisation)
    r[t] <- GetInterimResponse(selected_action, response_prob_a1_a, response_prob_a1_b)
    belief_a1 <- UpdateCurrentBelief(c(i,j,k,l), a1[t], r[t])
    i <- belief_a1[1]; j <- belief_a1[2]; k <- belief_a1[3]; l <- belief_a1[4]
  }

  O2_group_size <- belief_a1
  patient_index_list <- list(which(a1==0)[which(r[which(a1==0)]==1)], which(a1==0)[which(r[which(a1==0)]==0)],
                             which(a1==1)[which(r[which(a1==1)]==1)], which(a1==1)[which(r[which(a1==1)]==0)])

  # DTR Step 2 Inference
  for (group in 1:4) {
    group_size <- O2_group_size[group] - 1
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

  # Print summary
  cat("SIMULATION", sim, "\n")
  cat("DTR step 1 Summary: \n")
  cat("O1 |", bayes_prior_a1, "\n")
  cat("A1 |", "treatment a =", length(which(a1==0)), "; treatment b =", length(which(a1==1)), "\n")
  cat("O2 |", O2_group_size, "\n")
  cat("DTR step 2 Summary: \n")
  cat("A2 |", a2, "\n")
  cat("Y  |", y, "\n")
}