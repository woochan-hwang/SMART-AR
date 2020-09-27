# Wrapper for calling julia backend in main.R script

# LOAD DEPENDENCY ---------------------------
library(JuliaConnectoR)
juliaEval('include("/Users/WoochanH/Desktop/Projects/R_SMART-AR/code/Nested-CRDP/crdp-julia/nested_crdp.jl")' )

# LOAD HYPERPARAMETERS ----------------------
juliaLet("global N_PATIENTS = RtoJuliaVar", RtoJuliaVar = N_PATIENTS)
juliaLet("global bayes_prior_a1 = RtoJuliaVar", RtoJuliaVar = bayes_prior_a1)
juliaLet("global bayes_prior_a2 = RtoJuliaVar", RtoJuliaVar = bayes_prior_a2)
juliaLet("global degree_of_randomisation = RtoJuliaVar", RtoJuliaVar = degree_of_randomisation)
juliaLet("global degree_of_constraint = RtoJuliaVar", RtoJuliaVar = degree_of_constraint)

# FUNCTIONS ---------------------------------
JuliaIndex <- function(n,i,j,k) {
  juliaLet("global n = rtojuliavar", rtojuliavar = n)
  juliaLet("global i = rtojuliavar", rtojuliavar = i)
  juliaLet("global j = rtojuliavar", rtojuliavar = j)
  juliaLet("global k = rtojuliavar", rtojuliavar = k)
  index <- juliaEval("abs(DP_2_lin_index(n, i, j, k, 0))")
  # TODO needs bug fix for index returnng 0
  if (index == 0) {index <- 1}
  return (index)
}
# TODO adjust indexing difference for optimal action choice in julia vs R implementation
# TODO significant proportion on GC time during DTR step 2 julia version
# TODO refractoring of linear indexing
