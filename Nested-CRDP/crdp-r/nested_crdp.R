# Applying constrained randomised dynamic programming to SMART
# Woochan H. 2020, supervised by Dr David Robertson & Dr Sofia VIllar
# Includes functions for training and inference for simulations

# NOTES TO SELF
# All nested CRDPs should have equal degrees of randomisation and constraints given that ethically we should not be
# experimenting if we believed a certain arm to be definitively inferior to the others.

# Load dependencies ------------------------
source("crdp.R")
library(assertthat)

# KEY TRAINING FUNCTIONS -----------------------------
# Enumeration of possible state space for a binary trial given sample size
CreateValueMatrix <- function(sample_size, constraint) {

  assert_that(constraint <= 0.5)
  c <- as.integer(sample_size * constraint)

  value_matrix <- array(0, dim = c(sample_size+1, sample_size+1, sample_size+1))
  index_pairs <- c(NA, NA, NA, NA)  # return pairs that satisfy constraint
  t <- sample_size + 4

  for (i in 1:(t-3)){
    for (j in 1:(t-i-2)){
      for (k in 1:(t-i-j-1)){
        l <- (t-i-j-k)
        if (i+j < c) value_matrix[i,j,k] <- -sample_size
        if (k+l < c) value_matrix[i,j,k] <- -sample_size
        if (value_matrix[i,j,k] == 0) {
          index_pairs <- rbind(index_pairs, c(i,j,k,l))
        }
      }
    }
  }
  return (list("VALUE"=value_matrix, "INDICES"=index_pairs[-1,]))
}


# Create table with maximum IPTW output for each CRDP results of the second stage DTR
CreateMaxQ2Table <- function(q1_index_matrix, q2_list) {
  # create dataframe
  n1 <- dim(q1_index_matrix)[1]
  A2_max <- matrix(rep(0, (n1*4)), nrow=n1)
  A2_iptw <- matrix(rep(0, (n1*3)), nrow=n1)
  mydata <- data.frame(q1_index_matrix, A2_max, A2_iptw)
  names(mydata) <- c("A1_i", "A1_j", "A1_k", "A1_l", "A2_i_val", "A2_j_val", "A2_k_val", "A2_l_val",
                     "A2_IPTW_ij", "A2_IPTW_kl", "A2_IPTW_val")
  # fill with max_q2_val
  for (i in c("i", "j", "k", "l")) {
    column <- mydata[paste("A1", i, sep="_")]
    for (j in 1:n1) {
      n <- column[j,1]
      if (n >= 4) {
        val <- q2_list[[paste("N", n, "MAX", sep="_")]]
        ref <- paste("A2", i, "val", sep="_")
        mydata[ref][j,1] <- val
      }
    }
  }
  # IPTW
  mydata$"A2_IPTW_ij" <- (mydata$"A1_i"+mydata$"A1_j")*((mydata$"A2_i_val"/mydata$"A1_i")+(mydata$"A2_j_val"/mydata$"A1_j"))/2
  mydata$"A2_IPTW_kl" <- (mydata$"A1_k"+mydata$"A1_l")*((mydata$"A2_k_val"/mydata$"A1_k")+(mydata$"A2_l_val"/mydata$"A1_l"))/2
  # max_value
  larger_ij <- which(mydata$"A2_IPTW_ij" > mydata$"A2_IPTW_kl")
  larger_kl <- which(mydata$"A2_IPTW_ij" <= mydata$"A2_IPTW_kl")
  mydata$"A2_IPTW_val"[larger_ij] <- mydata$"A2_IPTW_ij"[larger_ij]
  mydata$"A2_IPTW_val"[larger_kl] <- mydata$"A2_IPTW_kl"[larger_kl]

  return (mydata)
}

# Backwards induction for Q1 matrix
BackwardsQ1Matrix <- function(sample_size, value_matrix, prior_matrix, p) {

  n      <- sample_size
  V      <- value_matrix
  Action <- array(NA, dim = c(n+1, n+1, n+1, n+1 ))
  prior.i <- prior_matrix[1]
  prior.j <- prior_matrix[2]
  prior.k <- prior_matrix[3]
  prior.l <- prior_matrix[4]

  for(t in (n+3):4){
    for (i in 1:(t-3)){
      for (j in 1:(t-i-2)){
        for (k in 1:(t-i-j-1)){

          l <- (t-i-j-k)

          expected.prob.c <- (i - 1 + prior.i) / (i - 1 + prior.i + j - 1 + prior.j)
          expected.prob.n <- (k - 1 + prior.k) / (k - 1 + prior.k + l - 1 + prior.l)

          Vcontrol <- expected.prob.c*(1 + V[i+1, j, k]) + (1 - expected.prob.c)*(0 + V[i, j+1, k])
          Vnovel   <- expected.prob.n*(1 + V[i, j, k+1]) + (1 - expected.prob.n)*(0 + V[i, j, k])

          if (p*Vcontrol + (1-p)*Vnovel > (1-p)*Vcontrol + p*Vnovel) {
            Action[n-(t-4), i, j, k] <- 0
            # Action 1 is optimal for the next patient.
          }
          if (p*Vcontrol + (1-p)*Vnovel < (1-p)*Vcontrol + p*Vnovel) {
            Action[n-(t-4), i, j, k] <- 1
            # Action 2 is optimal for the next patient.
          } else {
            if (p*Vcontrol + (1-p)*Vnovel == (1-p)*Vcontrol + p*Vnovel){
              Action[n-(t-4), i ,j, k] <- 2
              # Either Action 1 or 2 is optimal for the next patient.
            }
          }
          V[i,j,k] <- max( p*Vcontrol + (1-p)*Vnovel, (1-p)*Vcontrol + p*Vnovel )
        }
      }
    }
  }
  return(list("ACTION"=Action, "VALUE"=V))
}


# NESTED CRDP -------------------------------
# input: (sample_size, a1_prior_vector, a2_prior_vector, randomisation, constraint) where the prior vectors
# are vectors of length 4 (s_a,f_a,s_b,f_b); 0 <= constraint <= 0.5 <= randomisation <= 1; sample_size = int;
# return: a list with "Q1" and "Q2", where Q1 <- list("ACTION", "VALUE"); Q2 <- list("N_i_ACTION") for i in range(Q2)
TrainNestedCRDP <- function(sample_size, a1_prior_vector, a2_prior_vector, randomisation, constraint) {

  # DTR step 1 forward pass
  cat("[Running] Nested CRDP for DTR step 1 ... \n")
  ptm <- proc.time()
  q1 <- CreateValueMatrix(sample_size = sample_size, constraint = constraint)
  q1_value_matrix <- q1$VALUE
  q1_index_matrix <- q1$INDICES

  n1 <- dim(q1_index_matrix)[1]
  q1_max <- max(q1_index_matrix[,1])  # BUG: THIS NEEDS TO BE ACROSS ALL DIMS IN CASE OF NON-SYMMETRIC PRIORS
  q1_min <- max(min(q1_index_matrix[,1]), 4)
  cat("[Completed] Forward pass for DTR step 1:", (proc.time() - ptm)[[3]], "s \n")

  # DTR step 2 CRDP iteration
  q2_list <- list()
  cat("[Running] Nested CRDP for DTR step 2 ... \n")
  ptm <- proc.time()
  for (i in q1_min:q1_max) {
    q2 <- CRDP(n = i, prior_vector = a2_prior_vector, p = randomisation, Y = constraint)
    name1 <- paste("N", i, "ACTION", sep="_")
    name2 <- paste("N", i, "MAX", sep="_")
    q2_list[[name1]] <- q2$ACTION
    q2_list[[name2]] <- max(q2$VALUE)
  }
  cat("[Completed] Nested CRDP for DTR step 2:", (proc.time() - ptm)[[3]], "s \n")

  # create dataframe with maximum q2 val per q1 index
  cat("[Running] IPTW of step 2 values... \n")
  ptm <- proc.time()
  Q2Table <- CreateMaxQ2Table(q1_index_matrix, q2_list)
  cat("[Completed] IPTW value calculated for DTR step 2:", (proc.time() - ptm)[[3]], "s \n")
  cat("[Printing] Calculated Q2 value matrix sample... \n")
  print(head(Q2Table, 3))

  # Backwards induction
  cat("[Running] Backwards induction for step 1... \n")
  ptm <- proc.time()
  for (i in 1:n1) {
    r <- Q2Table[i,]
    q1_value_matrix[r[[1]],r[[2]],r[[3]]] <- r$"A2_IPTW_val"
  }
  q1_optim <- BackwardsQ1Matrix(sample_size, q1_value_matrix, a1_prior_vector, randomisation)
  cat("[Completed] Backwards induction for step 1:", (proc.time() - ptm)[[3]], "s \n")

  return (list("Q1"=q1_optim, "Q2"=q2_list))
}

