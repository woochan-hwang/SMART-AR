# DP optimisation code adopted from the binary bandit julia library

const BB_numerical_precision_64 = 1.0e-13
const BB_numerical_precision_32 = 1.0e-4

# COMMON FUNCTIONS -------------------------
# Converts a 4D state to linear index
function DP_2_lin_index( number_of_allocations , number_of_successes_arm_1 , number_of_failures_arm_1 , number_of_successes_arm_2 , number_of_remaining_allocations )
# Linear index equal to value_to_go_index used in policy loop when remaining_allocations = 1

    return div( number_of_allocations * ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - ( number_of_allocations - number_of_remaining_allocations + 1 ) * ( number_of_allocations - number_of_remaining_allocations + 2 ) * ( number_of_allocations - number_of_remaining_allocations + 3 ) * ( number_of_allocations - number_of_remaining_allocations + 4 ) , 24 ) + div( ( number_of_allocations - number_of_remaining_allocations + 1 ) * ( number_of_allocations - number_of_remaining_allocations + 2 ) * ( number_of_allocations - number_of_remaining_allocations + 3 ) - ( number_of_allocations - number_of_remaining_allocations - number_of_successes_arm_2 + 1 ) * ( number_of_allocations - number_of_remaining_allocations - number_of_successes_arm_2 + 2 ) * ( number_of_allocations - number_of_remaining_allocations - number_of_successes_arm_2 + 3 ) , 6 ) + div( ( number_of_allocations - number_of_remaining_allocations - number_of_successes_arm_2 + 1 ) * ( number_of_allocations - number_of_remaining_allocations - number_of_successes_arm_2 + 2 ) - ( number_of_allocations - number_of_remaining_allocations - number_of_successes_arm_2 - number_of_failures_arm_1 + 1 ) * ( number_of_allocations - number_of_remaining_allocations - number_of_successes_arm_2 - number_of_failures_arm_1 + 2 ) , 2 ) + number_of_successes_arm_1 + 1

end


# Converts hex encoded action sets of 4 to human readable format
function DP_policy_decoder( policy_bin :: Array{ UInt8 , 1 } )
# Each action can be recovered by running modulus of 4
# Minus 1 as optimal action encoding is different (action 0 in R == action 1 in Julia)
    decoded_policy = []
    for action_set in policy_bin
        hex = Int16(action_set)
        for i in 1:4
            hex = hex % 4
            push!(decoded_policy, hex - 1)
        end
    end

    return decoded_policy
end


# Apply constraint penalty and modify default value array
function ApplyConstraint!(value_to_go, number_of_allocations, constraint, prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2)

    c = constraint
    constraint_counter = [0,0]  # number allocated to [arm 1, arm 2] respectively
    number_of_observed_responses = number_of_allocations - 1
    value_to_go_lin_index = 0

    for number_of_successes_arm_2 = 0 : number_of_observed_responses , number_of_failures_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 ) , number_of_successes_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 - number_of_failures_arm_1 )
    @inbounds begin

        value_to_go_lin_index += 1
        constraint_counter[1] = number_of_successes_arm_1 + number_of_failures_arm_1 + prior_success_arm_1 + prior_failure_arm_1
        constraint_counter[2] = number_of_observed_responses - (number_of_successes_arm_1 + number_of_failures_arm_1) + prior_success_arm_2 + prior_failure_arm_2
        # applying penalty for non-constrained cases
        if ( constraint_counter[1] < c ) || ( constraint_counter[2] < c )
            value_to_go[ value_to_go_lin_index ] = - number_of_allocations
        end
    end #@inbounds
    end
end


# OFFLINE USAGE ------------------------------
# Output is the policy (i.e., actions for all states) with linear indices in binary encoding and the Bayes-expected number of successes
function CRDP_policy( number_of_allocations :: Int64 , degree_of_randomisation :: Float64, constraint :: Float64, prior_vector ) #

    prior_success_arm_1 = Int64(prior_vector[1])
    prior_failure_arm_1 = Int64(prior_vector[2])
    prior_success_arm_2 = Int64(prior_vector[3])
    prior_failure_arm_2 = Int64(prior_vector[4])

    value_to_go_32 :: Array{ Float32 , 1 } = zeros( Float32 , div( ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 , 6 ) + 1 )
    ApplyConstraint!(value_to_go_32, number_of_allocations, constraint, prior_success_arm_1 , prior_failure_arm_1, prior_success_arm_2, prior_failure_arm_2)
    action_matrix, value_matrix = CRDP_2_policy_bin_lin_with_finale( number_of_allocations, degree_of_randomisation, value_to_go_32, prior_success_arm_1, prior_failure_arm_1, prior_success_arm_2, prior_failure_arm_2)

    # returns decoded policy, value_matrix[1] is the expected bayes number of success following the optimal policy
    return DP_policy_decoder(action_matrix), value_matrix
end


function CRDP_2_policy_bin_lin_with_finale( number_of_allocations :: Int64 , degree_of_randomisation :: Float64, value_to_go :: Array{ Float32 , 1 } , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )

    # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations ) not needed as it is given by value_to_go
    p = degree_of_randomisation
    exp_value_of_action_1 = 0.0 # needed after the for loop
    exp_value_of_action_2 = 0.0 # needed after the for loop

    # backwards recursion: t-th step
    action_bin = zeros( UInt8 , div( number_of_allocations * ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 * 4 , 96 ) + 1 ) # this is the correct size of the action in binary encoding
    mod_lin = mod( div( number_of_allocations * ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 * 4 , 24 ) + 1 , 4 ) # modulo of size of action_lin divided by 4
    action_bin_temp :: UInt8 = 0x00 # encoded up to 4 actions to be written at the next position of action_bin
    action_bin_temp_index = mod( 4 - mod_lin , 4 ) # number of actions currently encoded in action_bin_temp; initialise (up to 3 first values may be unused) so that the final value in action_lib refers to the initial state
    lin_index = 0 # last filled position of vector action_bin
    for number_of_observed_responses = ( number_of_allocations - 1 ) : -1 : 0
        value_to_go_lin_index = 0
        for number_of_successes_arm_2 = 0 : number_of_observed_responses , number_of_failures_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 ) , number_of_successes_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 - number_of_failures_arm_1 )
        @inbounds begin

            value_to_go_lin_index += 1

            # current belief of each arm given prior
            belief_of_success_arm_1 = ( ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 ) )
            belief_of_success_arm_2 = ( ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) )

            # value given DP method
            value_to_go_if_arm_1 = belief_of_success_arm_1 * ( 1.0 + value_to_go[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) + number_of_failures_arm_1 + 1 ] ) + ( 1.0 - belief_of_success_arm_1 ) * value_to_go[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) + ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) ]
            value_to_go_if_arm_2 = belief_of_success_arm_2 * ( 1.0 + value_to_go[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) , 2 ) ] ) + ( 1.0 - belief_of_success_arm_2 ) * value_to_go[ value_to_go_lin_index + number_of_failures_arm_1 + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) ]

            # expectation of value given degree of randomisation
            exp_value_of_action_1 = ( p * value_to_go_if_arm_1 ) + ( ( 1 - p ) * value_to_go_if_arm_2 )
            exp_value_of_action_2 = ( ( 1 - p ) * value_to_go_if_arm_1 ) + ( p * value_to_go_if_arm_2 )

            # update value_to_go array
            if ( exp_value_of_action_1 - exp_value_of_action_2 ) > BB_numerical_precision_32 * ( exp_value_of_action_1 + exp_value_of_action_2 )
                value_to_go[ value_to_go_lin_index ] += Float32( exp_value_of_action_1 )
                action_bin_temp = ( action_bin_temp * 0x04 ) + 0x01 # action 1 is optimal, action 2 is not optimal
            elseif ( exp_value_of_action_2 - exp_value_of_action_1 ) > BB_numerical_precision_32 * ( exp_value_of_action_1 + exp_value_of_action_2 )
                value_to_go[ value_to_go_lin_index ] += Float32( exp_value_of_action_2 )
                action_bin_temp = ( action_bin_temp * 0x04 ) + 0x02 # action 1 is not optimal, action 2 is optimal
            else #if exp_value_of_action_1 approx== exp_value_of_action_2
                value_to_go[ value_to_go_lin_index ] += Float32( ( exp_value_of_action_1 + exp_value_of_action_2 ) / 2 )
                action_bin_temp = ( action_bin_temp * 0x04 ) + 0x03 # randomise between actions 1 and 2 as both actions are optimal
            end

            # save action
            action_bin_temp_index += 1
            if action_bin_temp_index == 4
                lin_index += 1
                action_bin[ lin_index ] = action_bin_temp
                action_bin_temp = 0x00
                action_bin_temp_index = 0
            end

        end # @inbounds
        end
    end
    # complete action matrix, value_to_go[1] is the expected bayes number of success following the optimal policy
    return action_bin , value_to_go
end


# ONLINE USAGE ------------------------------
# Output is the immediate action and the Bayes-expected number of successes (as Float64 irrespectively of the value of "precision")
# Uses linear indexing
function CRDP_online( number_of_allocations :: Int64 , degree_of_randomisation :: Float64, constraint :: Float64, prior_vector, float_version :: Int64 = Int64( 32 ))

    prior_success_arm_1 = Int64(prior_vector[1])
    prior_failure_arm_1 = Int64(prior_vector[2])
    prior_success_arm_2 = Int64(prior_vector[3])
    prior_failure_arm_2 = Int64(prior_vector[4])

    if float_version == 32
        # Create value array and apply constraint
        value_to_go_32 :: Array{ Float32 , 1 } = zeros( Float32 , div( ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 , 6 ) + 1 )
        ApplyConstraint!(value_to_go_32, number_of_allocations, constraint, prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2)
        # Run CRDP on line algorithm
        CRDP_2_action_lin_with_finale( number_of_allocations , degree_of_randomisation, value_to_go_32 , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )
    else
        println("undefined")
        #TODO : Implement float16 and float64 version
    end
end

# Float32 version of value_to_go
function CRDP_2_action_lin_with_finale( number_of_allocations :: Int64 , degree_of_randomisation :: Float64, value_to_go :: Array{ Float32 , 1 } , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )

    # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations ) not needed as it is given by value_to_go
    p = degree_of_randomisation
    # backwards recursion: t-th step
    exp_value_of_action_1 = 0.0 # needed after the for loop
    exp_value_of_action_2 = 0.0 # needed after the for loop

    for number_of_observed_responses = ( number_of_allocations - 1 ) : -1 : 0
        value_to_go_lin_index = 0
        for number_of_successes_arm_2 = 0 : number_of_observed_responses , number_of_failures_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 ) , number_of_successes_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 - number_of_failures_arm_1 )
        @inbounds begin

            value_to_go_lin_index += 1

            # current belief of each arm given prior
            belief_of_success_arm_1 = ( ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 ) )
            belief_of_success_arm_2 = ( ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) )

            # value given DP method
            value_to_go_if_arm_1 = belief_of_success_arm_1 * ( 1.0 + value_to_go[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) + number_of_failures_arm_1 + 1 ] ) + ( 1.0 - belief_of_success_arm_1 ) * value_to_go[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) + ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) ]
            value_to_go_if_arm_2 = belief_of_success_arm_2 * ( 1.0 + value_to_go[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) , 2 ) ] ) + ( 1.0 - belief_of_success_arm_2 ) * value_to_go[ value_to_go_lin_index + number_of_failures_arm_1 + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) ]

            # expectation of value given degree of randomisation
            exp_value_of_action_1 = ( p * value_to_go_if_arm_1 ) + ( ( 1 - p ) * value_to_go_if_arm_2 )
            exp_value_of_action_2 = ( ( 1 - p ) * value_to_go_if_arm_1 ) + ( p * value_to_go_if_arm_2 )

            # update value_to_go array
            if ( exp_value_of_action_1 - exp_value_of_action_2 ) > BB_numerical_precision_32 * ( exp_value_of_action_1 + exp_value_of_action_2 )
                value_to_go[ value_to_go_lin_index ] += Float32( exp_value_of_action_1 )
            elseif ( exp_value_of_action_2 - exp_value_of_action_1 ) > BB_numerical_precision_32 * ( exp_value_of_action_1 + exp_value_of_action_2 )
                value_to_go[ value_to_go_lin_index ] += Float32( exp_value_of_action_2 )
            else #if exp_value_of_action_1 approx== exp_value_of_action_2
                value_to_go[ value_to_go_lin_index ] += Float32( ( exp_value_of_action_1 + exp_value_of_action_2 ) / 2 )
            end

        end # @inbounds
        end
    end

    if ( exp_value_of_action_1 - exp_value_of_action_2 ) > BB_numerical_precision_32 * ( exp_value_of_action_1 + exp_value_of_action_2 )
        action = Int8( 1 ) # action 1 (i.e., arm 1 with probability of p)
    elseif ( exp_value_of_action_2 - exp_value_of_action_1 ) > BB_numerical_precision_32 * ( exp_value_of_action_1 + exp_value_of_action_2 )
        action = Int8( 2 ) # action 2 (i.e., arm 2 with probability of p)
    else
        action = Int8( 3 ) # randomise between actions 1 and 2
    end

    return action , Float64(value_to_go[ 1 ])
end

# Example of running simulations using online version of CRDP
function OnlineCRDP(N, degree_of_randomisation, degree_of_constraint, prior_vector)

    constraint = degree_of_constraint * N

    for number_allocated in [N:-1:1;]

        println("RESULTS for $number_allocated")

        horizon = number_allocated
        action, value = CRDP_online(horizon , degree_of_randomisation, constraint, prior_vector)

        println("Action: $action for value $value")

        treatment = SelectTreatment(action, degree_of_randomisation)  # randomised DP
        response = TreatmentResponse(treatment, success_prob_action_1, success_prob_action_2)  # get response from Bernoulli distribution

        println("treatment: $treatment, response: $response")

        UpdateCurrentBelief!(prior_vector, treatment, response)  # update prior vector

        println("new prior vector: $prior_vector")
    end

end
