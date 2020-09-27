# Nested CRDP
  
Woochan H. 2020, supervised by Dr David Robertson & Dr Sofia VIllar

### SUMMARY 

This section provides code for applying constrainted randomised dynamic programming (CRDP) [1] in the context of a two stage dynamic treatment regime. Adaptive randomisation [2] has been explored for a long time in the context of single stage randomised control trials but few attempts have been made to assess the effects for sequential study designs. 

In this project we extend the CRDP algorithm to a two-stage, two-armed scenario (Nested CRDP) and will be assessing the merits of the design using type 1 error, statistical power and relative patient benefit. This will be compared to the SMART-AR algorithm [3] proposed by Cheung et al. using the same metrics. One advantage of our approach is the simplicity of tuning the parameters. 

The original R code for a single stage CRDP provided by Faye et al. is in the 'crdp.R' script. 

The extension of the constrained randomised dynamic programming method to a two stage binary treatment setting is contained in 'nested_crdp.R'. The training and simulation is run from 'main.R' and all experiment parameters are maintained at the top of this script. Most of the notations are inherited from the original crdp code to increase readability. We provide a julia implementation which showed significant improvement in performance (memory and computation speed), allowing for larger simulations. The backend can be specified as a hyperparameter in the "main.R" script. 

The current version is able to run up to a sample size of (75 using the R version; 250 using the julia version) on a standard laptop (~5min for enumeration & 10,000 simulations). The setting for comparing the first stage DTR regardless of the following outcome ( null_hypothesis <- “a1”) is equivalent to the original 1 stage CRDP setting. The constraint and randomisation parameters follow the recommendation from the original paper but will need to be further validated for optimality in the DTR scenario. 


**Notes to self
1. SMART-AR calibrates based on average BDI improvement -> or do we want to maximise the proportion allocated to superior treatment


### TODO
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

#### References
 
[1] S. Faye, et al. "A Bayesian adaptive design for clinical trials in rare diseases." Computational statistics & data analysis 113 (2017): 136-153.

[2] Robertson, David S., et al. "Response-adaptive randomization in clinical trials: from myths to practical considerations." arXiv preprint arXiv:2005.00564 (2020).

[3] Cheung, Ying Kuen, Bibhas Chakraborty, and Karina W. Davidson. "Sequential multiple assignment randomized trial (SMART) with adaptive randomization for quality improvement in depression treatment program." Biometrics 71.2 (2015): 450-459.