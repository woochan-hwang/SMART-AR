# Applying constrained randomised dynamic programming to SMART julia version
# Woochan H. 2020, supervised by Dr David Robertson & Dr Sofia VIllar

# LOAD DEPENDENCY ---------------------------
using Distributions
include("utility_functions.jl")
include("crdp.jl")


# KEY TRAINING FUNCTIONS --------------------
# Enumeration of possible state space for a binary trial given sample size
function CreateValueMatrix(sample_size, constraint)

    @assert constraint <= 0.5
    c = Int64(floor(sample_size * constraint))

    value_matrix :: Array{ Float32 , 1 } = zeros( Float32 , div( ( sample_size + 1 ) * ( sample_size + 2 ) * ( sample_size + 3 ) - 1 * 2 * 3 , 6 ) + 1 )
    value_matrix_lin_index = 0
    index_pairs = []
    max_allocated = 0

    for successes_arm_2 = 0 : sample_size, failures_arm_1 = 0 : ( sample_size - successes_arm_2) , successes_arm_1 = 0 : ( sample_size - successes_arm_2 - failures_arm_1)
    @inbounds begin

        value_matrix_lin_index += 1
        if c <= (successes_arm_1 + failures_arm_1) <= sample_size - c
            # store index as constraint satisfied
            push!(index_pairs, [successes_arm_1, failures_arm_1, successes_arm_2])
        else
            # penalty as constraint not satisfied
            value_matrix[ value_matrix_lin_index ] = - sample_size
        end
    end #@inbounds
    end

    return value_matrix, index_pairs
end


# Create table with maximum IPTW output for each CRDP results of the second stage DTR
function CreateMaxQ2Table(sample_size, q1_index_matrix, q2_list)

    n1 = size(q1_index_matrix)[1]
    mydata = zeros( Float32, n1, 11)
    mydata[:,1:3] = reshape(collect(Iterators.flatten(q1_index_matrix)), (n1, 3))
    mydata[:,4] = sample_size .- (mydata[:,1] + mydata[:,2] + mydata[:,3])
    for column = 1:4, row = 1:n1
        n = max(Int64(mydata[row, column]), 4) # minimum horizon bounded to 4
        if n > 86
        end
        mydata[row, column + 4] = q2_list[n][1] # q2 value when horizon is n
    end

    # IPTW
    mydata[:, 9] = (mydata[:,1] + mydata[:,2]) .* ((mydata[:,5] ./ mydata[:,1]) + (mydata[:,6] ./ mydata[:,2])) ./ 2
    mydata[:, 10] = (mydata[:,3] + mydata[:,4]) .* ((mydata[:,7] ./ mydata[:,3]) + (mydata[:,8] ./ mydata[:,4])) ./ 2
    # max_value
    larger_ij = findall(mydata[:,9] .> mydata[:,10])
    larger_kl = findall(mydata[:,9] .< mydata[:,10])
    mydata[larger_ij, 11] = mydata[larger_ij, 9]
    mydata[larger_kl, 11] = mydata[larger_kl, 10]

    return mydata
end

# Backwards induction for Q1 matrix
function BackwardsQ1Matrix(sample_size, value_matrix, prior_vector, randomisation)

    prior_success_arm_1 = Int64(prior_vector[1])
    prior_failure_arm_1 = Int64(prior_vector[2])
    prior_success_arm_2 = Int64(prior_vector[3])
    prior_failure_arm_2 = Int64(prior_vector[4])

    # Run CRDP with updated value matrix based on q2 findings
    action_matrix, value_matrix = CRDP_2_policy_bin_lin_with_finale(sample_size, randomisation, value_matrix, prior_success_arm_1, prior_failure_arm_1, prior_success_arm_2, prior_failure_arm_2)
    # return policy for Q1
    return Dict("ACTION" => DP_policy_decoder(action_matrix), "VALUE" => value_matrix)
end


# NESTED CRDP -------------------------------
# input: (sample_size, a1_prior_vector, a2_prior_vector, randomisation, constraint) where the prior vectors
# are vectors of length 4 (s_a,f_a,s_b,f_b); 0 <= constraint <= 0.5 <= randomisation <= 1; sample_size = int;
# return: a list with "Q1" and "Q2", where Q1 <- list("ACTION", "VALUE"); Q2 <- list("N_i_ACTION") for i in range(Q2)
function TrainNestedCRDP(sample_size, a1_prior_vector, a2_prior_vector, randomisation, constraint)

    sample_size = Int64(sample_size)
    randomisation = Float64(randomisation)
    constraint = Float64(constraint)

    # DTR step 1 forward pass
    println("[Running] Nested CRDP for DTR step 1 ...")
    q1_value_matrix, q1_index_matrix = @time CreateValueMatrix(sample_size, constraint)
    q1_max = max(maximum(collect(Iterators.flatten(q1_index_matrix))), sample_size - minimum(collect(Iterators.flatten(q1_index_matrix))))
    # TODO check q1_min ?-1
    println("q1_max: ", q1_max)
    println("[Completed] Forward pass for DTR step 1")

    # DTR step 2 CRDP iteration
    q2_action_dict = Dict() # dict as policy is returned and called in train.R script
    q2_value_list = [] # q2_list[index] where index == q2 horizon
    println("[Running] Nested CRDP for DTR step 2 ... ")
    # TODO check indexing for q1_max and q1_index_array
    for horizon in 1 : q1_max + 1
        print("horizon: ", horizon)
        # TODO change to offline policy bin version and update datastructure for q2_list
        action, value = @time CRDP_policy(horizon, randomisation, constraint, a2_prior_vector)
        push!(q2_value_list, value)
        q2_action_dict[string("N_", horizon, "_ACTION")] = action
    end
    println("[Completed] Nested CRDP for DTR step 2")

    # Create dataframe with maximum q2 value per q1 index
    println("[Running] IPTW of step 2 values...")
    q2_table = @time CreateMaxQ2Table(sample_size, q1_index_matrix, q2_value_list)
    println("[Completed] IPTW value calculated for DTR step 2")

    # Backwards induction
    indices = []
    println("[Running] Backwards induction for step 1...")
    for i in 1 : size(q1_index_matrix)[1]
        r = q2_table[i,:]
        lin_index = Int64(DP_2_lin_index(sample_size, r[1], r[2], r[3], 1))
        q1_value_matrix[lin_index] = r[11]
    end
    q2_table = nothing
    q1_optim = @time BackwardsQ1Matrix(sample_size, q1_value_matrix, a1_prior_vector, randomisation)
    println("[Completed] Backwards induction for step 1")

    return Dict("Q1" => q1_optim, "Q2" => q2_action_dict, "Q2_Keys" => keys(q2_action_dict))
end
