# Nested CRDP

The original code provided by Faye et al. is in the 'crdp.R' script. 

S. Faye, et al. "A Bayesian adaptive design for clinical trials in rare diseases." Computational statistics & data analysis 113 (2017): 136-153.

The extension of the constrained randomised dynamic programming method to a 2-stage binary dynamic treatment regime setting is contained in 'nested_crdp.R'. The training and simulation is run from 'train.R' and all experiment parameters are maintained at the top of this script. Most of the notations are inherited from the original crdp code to increase readability. 

The current version is able to run up to a sample size of 75 on a standard laptop (~5min for enumeration & 10,000 simulations). The setting for comparing the first stage DTR regardless of the following outcome ( null_hypothesis <- “a1”) is equivalent to the original 1 stage CRDP setting. The constraint and randomisation parameters follow the recommendation from the original paper but will need to be further validated for optimality in the DTR scenario. 