# Applying constrained randomised dynamic programming to SMART
# Woochan H. 2020, supervised by Dr David Robertson & Dr Sofia VIllar
# Includes misc functions for hypothesis testing of simulation outputs


# LOAD DEPENDENCY ---------------------------
using Distributions
using HypothesisTests
include("crdp.jl")

# KEY INFERENCE FUNCTIONS ---------------------------
function BernoulliResponse(probability)

     d = Distributions.Binomial(1, probability)
     return rand(d)

end

function GetOptimalAction(optimal_policy, total_allocated, n, i, j, k)

    index = Int64(DP_2_lin_index(total_allocated, i, j, k, 1))
    #println("INDEX: ", index, " where: [total, n, i, j, k] = [$total_allocated, $n, $i, $j, $k]")
    return optimal_policy[index]

end

function SelectTreatment(optimal_action, degree_of_randomisation)

    if optimal_action == 0
        selected_treatment = 1 - BernoulliResponse(degree_of_randomisation)
    elseif optimal_action == 1
        selected_treatment = BernoulliResponse(degree_of_randomisation)
    else
        selected_treatment = BernoulliResponse(0.5)
    end
    return selected_treatment

end

function GetInterimResponse(action, success_prob_action_1, success_prob_action_2)

    if action == 0
        response = BernoulliResponse(success_prob_action_1)
    elseif action == 1
        response = BernoulliResponse(success_prob_action_2)
    else
        response = BernoulliResponse(0.5)
    end
    return response

end

function UpdateCurrentBelief!(prior_vector, treatment, response)

    if treatment == 0 && response == 1
        prior_vector[1] += 1
    elseif treatment == 0 && response == 0
        prior_vector[2] += 1
    elseif treatment == 1 && response == 1
        prior_vector[3] += 1
    elseif treatment == 1 && response == 0
        prior_vector[4] += 1
    end

end

function GetOutcome(a1, a2, response_prob_a2_aa, response_prob_a2_ab, response_prob_a2_ba, response_prob_a2_bb)

    if a2 == 0
        if a1 == 0
            response = BernoulliResponse(response_prob_a2_aa)
        elseif a2 == 1
            response = BernoulliResponse(response_prob_a2_ba)
        end
    elseif a2 == 1
        if a1 == 0
            response = BernoulliResponse(response_prob_a2_ab)
        elseif a2 == 1
            response = BernoulliResponse(response_prob_a2_bb)
        end
    end
    return response

end



# PERFORMANCE METRICS -----------------------
# Fisher's Exact Test, small samples of two binomial distributions
function HypothesisTesting(dataframe, null, alternative, a1_h = nothing, r_h = nothing)

  null == "a2.given.a1" && @assert !isnothing(a1_h)
  null == "a2.given.r" && @assert (!isnothing(a1_h) && !isnothing(r_h))
  a1 = dataframe["A1"]
  r = dataframe["R"]
  a2 = dataframe["A2"]
  y = dataframe["Y"]

  if null == "a1"
    sim_results = [length(findall(r[a1.==0].==1)), length(findall(r[a1.==0].==0)), length(findall(r[a1.==1].==1)), length(findall(r[a1.==1].==0))]
  elseif null == "a2.combined"
    sim_results = [length(findall(y[a2.==0].==1)), length(findall(y[a2.==0].==0)), length(findall(y[a2.==1].==1)), length(findall(y[a2.==1].==0))]
  elseif null == "a2.given.a1"
    sim_results = [length(findall(y[findall(a2[a1.==a1_h].==0)].==1)), length(findall(y[findall(a2[a1.==a1_h].==0)].==0)),
                         length(findall(y[findall(a2[a1.==a1_h].==1)].==1)), length(findall(y[findall(a2[a1.==a1_h].==1)].==0))]
  elseif null == "a2.given.r"
  end
  if sum(sim_results[1:2]) == 0 || sum(sim_results[3:4]) == 0 || sim_results[1] + sim_results[3] == 0 || sim_results[1] + sim_results[3] == 0
      # sample not sufficent for fisher's exact testing
      println("Null hypothesis rejected due to insufficient table for Fisher's Exact Test: ", sim_results)
      return 0.001
  else
      results = HypothesisTests.FisherExactTest(sim_results[1], sim_results[3], sim_results[2], sim_results[4])
      return HypothesisTests.pvalue(results)
  end
end

function PatientOutcome(dataframe, null, a1_h, r_h, superior_treatment)

  null == "a2.given.a1" && @assert !isnothing(a1_h)
  a1 = dataframe["A1"]
  a2 = dataframe["A2"]

  if null == "a1"
    outcome = length(findall(a1.==superior_treatment)) / length(a1)
  elseif null == "a2.combined"
    outcome = length(findall(a2.==superior_treatment)) / length(a2)
  elseif null == "a2.given.a1"
    outcome = length(findall(a2[a1.==a1_h].==superior_treatment)) / length(findall(a1.==a1_h))
  elseif null == "a2.given.r"
  end
  return outcome

end

function ExperimentSetting(null, response_prob_a1_a, response_prob_a1_b, response_prob_a2_aa, response_prob_a2_ab, response_prob_a2_ba, response_prob_a2_bb, a1_h)

  if a1_h == 0
      a1_h = 'a'
  else
      a1_h = 'b'
  end

  if null == "a1"
    prob1 = response_prob_a1_a
    prob2 = response_prob_a1_b
  elseif null == "a2.combined"
    prob1 = (response_prob_a2_aa + response_prob_a2_ba) / 2
    prob2 = (response_prob_a2_ab + response_prob_a2_bb) / 2
  elseif null == "a2.given.a1"
    prob1 = paste("response_prob_a2_", a1_h, "a", sep="")
    prob2 = paste("response_prob_a2_", a1_h, "b", sep="")
  elseif null == "a2.given.r"
    # look up interaction design
  end

  if prob1 == prob2
      ground_truth = true
  else
      ground_truth = false
  end

  if prob1 > prob2
      superior_treatment = 0
  else
      superior_treatment = 1
  end

  return Dict("ground_truth" => ground_truth, "superior_treatment" => superior_treatment)
end


# OUTPUT FUNCTIONS --------------------------

function PrintHyperParameters()

  listOfParameters = Dict()
  parameterNames = ["degree_of_randomisation", "degree_of_constraint", "bayes_prior_a1", "bayes_prior_a2",
                    "N_SIM", "N_PATIENTS", "response_prob_a1_a", "response_prob_a1_b", "response_prob_a2_aa",
                    "response_prob_a2_ab", "response_prob_a2_ba", "response_prob_a2_bb", "null_hypothesis",
                    "alternative", "a1_h", "r_h", "significance_level", "ground_truth", "superior_treatment"]
  for obj in parameterNames
    listOfParameters[obj] = eval(Meta.parse(obj))
    println(obj, ": ", listOfParameters[obj])
  end

  return listOfParameters

end