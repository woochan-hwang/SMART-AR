# Applying constrained randomised dynamic programming to a two stage binary outcome design
# Woochan H. 2020, supervised by Dr David Robertson & Dr Sofia Villar
# Julia optimised version of main.R, supports limited functionality at the moment.

# LOAD DEPENDENCY ---------------------------
using Random
include("nested_crdp.jl")
include("utility_functions.jl")


# HYPERPARAMETERS ---------------------------
# RAR algorithm design parameters
const degree_of_randomisation = 0.9  # DoR: prob of selecting arm with higher expectation; 0.5 <= DoR <= 1
const degree_of_constraint = 0.15  # DoC * N_patients = minimum allocation per DTR arm; DoC <= 0.5

# Prior probability; c(success_arm_a, failure_arm_a, success_arm_b, failure_arm_b)
const bayes_prior_a1 = [1,1,1,1]  # stage 1 prior
const bayes_prior_a2 = [1,1,1,1]  # stage 2 prior

# Simulation parameters
const N_SIM = 1000  # number of simulation repeats
const N_PATIENTS = 10

# True distribution parameters
const response_prob_a1_a = 0.8
const response_prob_a1_b = 0.5
const response_prob_a2_aa = 0.5  # P( Response == 1 | a1 = a, a2 = a )
const response_prob_a2_ab = 0.5  # P( Response == 1 | a1 = a, a2 = b )
const response_prob_a2_ba = 0.5
const response_prob_a2_bb = 0.5

# Hypothesis testing parameters
const null_hypothesis = "a1"  # c("a1", "a2.combined", "a2.given.a1", "a2.given.r")
const alternative = "both"  # c("both", "one")
const a1_h = 0; const r_h = 0  # required for "a2.given.a1", "a2.given.r"
const significance_level = 0.1  # set as 0.1 in original paper

# True outcome based on given parameters
ground_truth = ExperimentSetting(null_hypothesis, response_prob_a1_a, response_prob_a1_b, response_prob_a2_aa,
                                  response_prob_a2_ab, response_prob_a2_ba, response_prob_a2_bb, a1_h)["ground_truth"]
superior_treatment = ExperimentSetting(null_hypothesis, response_prob_a1_a, response_prob_a1_b, response_prob_a2_aa,
                                        response_prob_a2_ab, response_prob_a2_ba, response_prob_a2_bb, a1_h)["superior_treatment"]
parameterList = PrintHyperParameters()


# TRAIN OPTIMAL POLICY -------------------------
optimal_policy = TrainNestedCRDP(N_PATIENTS, bayes_prior_a1, bayes_prior_a2, degree_of_randomisation, degree_of_constraint)
optimal_q1_policy = optimal_policy["Q1"]["ACTION"]
optimal_q2_policy = optimal_policy["Q2"]

# RUN SIMULATIONS ------------------------------
Random.seed!(1234)
pval = Array{Float32}(undef, N_SIM)
patient_outcome = Array{Float32}(undef, N_SIM)
println("[Running] ", N_SIM, "simulation iterations for ", N_PATIENTS, " samples")

for sim = 1 : N_SIM

    belief_a1 = [0, 0, 0, 0]
    a1 = r = a2 = y = Array{Float32}(undef, N_PATIENTS)

    # DTR Step 1 Inference
    for t = 1 : N_PATIENTS
        n = N_PATIENTS - t + 1
        optimal_action = GetOptimalAction(optimal_q1_policy, N_PATIENTS, n, belief_a1[1], belief_a1[2], belief_a1[3])
        a1[t] = selected_action = SelectTreatment(optimal_action, degree_of_randomisation)
        r[t] = GetInterimResponse(selected_action, response_prob_a1_a, response_prob_a1_b)
        UpdateCurrentBelief!(belief_a1, a1[t], r[t])
    end

    patient_index_list = [findall(a1.==0)[findall(r[a1.==0].==1)], findall(a1.==0)[findall(r[a1.==0].==0)],
                          findall(a1.==1)[findall(r[a1.==1].==1)], findall(a1.==1)[findall(r[a1.==1].==0)]]

    # DTR Step 2 Inference
    for group = 1 : 4
        group_size = belief_a1[group]
        group_name = string("N_", group_size, "_ACTION")
        if group_size == 0
            continue  # skip iteration for groups with no patients
        elseif 0 < group_size < 4
            group_name = "N_4_ACTION"
        end
        optimal_group_policy = optimal_q2_policy[group_name]
        belief_a2 = [0, 0, 0, 0]

        # A2 selection based on optimal policy per subgroup
        for t in 1 : group_size
            patient_index = patient_index_list[group][t]
            n = group_size - t + 1
            optimal_action = GetOptimalAction(optimal_group_policy, group_size, n, belief_a2[1], belief_a2[2], belief_a2[3])
            a2[patient_index] = selected_action = SelectTreatment(optimal_action, degree_of_randomisation)
            y[patient_index] = GetOutcome(a1[patient_index], selected_action, response_prob_a2_aa, response_prob_a2_ab,
                                         response_prob_a2_ba, response_prob_a2_bb)
            UpdateCurrentBelief!(belief_a2, a2[patient_index], y[patient_index])
        end
    end

    df = Dict("A1" => a1, "R" => r, "A2" => a2, "Y" => y)
    pval[sim] = HypothesisTesting(df, null_hypothesis, alternative, a1_h, r_h)
    patient_outcome[sim] = PatientOutcome(df, null_hypothesis, a1_h, r_h, superior_treatment)
end
println("[Completed] ", N_SIM, " simulation iterations")


# HYPOTHESIS TESTING ------------------------
if ground_truth
  type_1_error = length(findall(pval .< significance_level)) / N_SIM
  println("Type 1 Error: ", type_1_error)
else
  power = length(findall(pval .< significance_level)) / N_SIM
  println("Statistical Power: ", power)
end

println("Proportion of patients allocated to superior arm: ", mean(patient_outcome))
