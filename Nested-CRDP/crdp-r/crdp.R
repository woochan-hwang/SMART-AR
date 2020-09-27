# Modified from code provided by Williamson, S. Faye, et al. "A Bayesian adaptive design for clinical trials
# in rare diseases." Computational statistics & data analysis 113 (2017): 136-153.

# NOTATION
# n: number of patients in trial.
# prior.i: prior number of successes on treatment A.
# prior.j: prior number of failures on treatment A.
# prior.k: prior number of successes on treatment B.
# prior.l: prior number of failures on treatment B.
# V: value function representing the maximum expected total reward (i.e. number of successes) in the rest of the trial after t patients have been treated.
# t: number of patients that have been treated.
# p: the degree of randomisation.
# Y: minimum number of observations required on each treatment arm (i.e. the degree of constraining).
# i: observed number of successes on treatment A.
# j: observed number of failures on treatment A.
# k: observed number of successes on treatment B.
# l: observed number of failures on treatment B.

CRDP <- function(n, prior_vector, p, Y){

  V      <- array(0, dim = c(n+1, n+1, n+1))
  Action <- array(0, dim = c(n+1, n+1, n+1, n+1 ))
  t      <- n+4
  prior.i <- prior_vector[1]
  prior.j <- prior_vector[2]
  prior.k <- prior_vector[3]
  prior.l <- prior_vector[4]
  Y <- as.integer(Y*n)

  # Create value matrix and assess constraint
  for (i in 1:(t-3)){
    for (j in 1:(t-i-2)){
      for (k in 1:(t-i-j-1)){

        l <- (t-i-j-k)
        if (i+j < Y) V[i,j,k] <- -n
        if (k+l < Y) V[i,j,k] <- -n
      }
    }
  }
  # Backwards induction
  for(t in (n+3):4){
    for (i in 1:(t-3)){
      for (j in 1:(t-i-2)){
        for (k in 1:(t-i-j-1)){

          l <- (t-i-j-k)

          expected.prob.c <- (i - 1 + prior.i) / (i - 1 + prior.i + j - 1 + prior.j)  # doesn't work for 0 priors
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
