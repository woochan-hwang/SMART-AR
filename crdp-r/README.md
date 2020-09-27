# Nested CRDP

The original code provided by Faye et al. is in the 'crdp.R' script. 

S. Faye, et al. "A Bayesian adaptive design for clinical trials in rare diseases." Computational statistics & data analysis 113 (2017): 136-153.

The extension of the constrained randomised dynamic programming method to a 2-stage binary dynamic treatment regime setting is contained in 'nested_crdp.R'. The training and simulation is run from 'train.R' and all experiment parameters are maintained at the top of this script. Most of the notations are inherited from the original crdp code to increase readability. 

The current version is able to run up to a sample size of 75 on a standard laptop (~5min for enumeration & 10,000 simulations). The setting for comparing the first stage DTR regardless of the following outcome ( null_hypothesis <- “a1”) is equivalent to the original 1 stage CRDP setting. The constraint and randomisation parameters follow the recommendation from the original paper but will need to be further validated for optimality in the DTR scenario. 


1. SMART-AR calibrates based on average BDI improvement -> or do we want to maximise the proportion allocated to superior treatment


TASKS
1. Design interactions between a1 and a2 -> likely require larger N, check how much memory is required. 
2. Discuss optimal DP, JULIA code with Sofia 
3. Include drop out rate as parameter and link it with patient benefit with different strengths. Assumption here is that 
participants who see a benefit are more likely to stay on the trial - especially in mHealth. 
 - This will have direct implications on study power
4. Consolidate what limitations Bayesian perspective may have on statistical power analysis of DTR. 
 - What makes a good method? 
 - SMART-AR method is currently very limited and calibrated to the available prior. 
 - Effect of priors on statistical power 
 - How much calibration do we really need? Is it fair / unfair to use a uncalibrated version? given practicality of use
 - SMART-AR currently performs worse on patient benefit 
 - Critique there method of proportion of finding optimal DTR vs statistical power
 
