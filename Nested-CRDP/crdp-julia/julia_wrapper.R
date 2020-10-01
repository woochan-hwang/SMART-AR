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
