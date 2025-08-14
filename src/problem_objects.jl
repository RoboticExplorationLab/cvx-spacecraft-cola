# Represents a satellite conjunction event between two space objects in the bplane (encounter frame).
struct ConjunctionEvent
    sat1_tca          # e.g. length-6 vector [r; v] at TCA
    sat2_tca          # e.g. length-6 vector [r; v] at TCA
    cov_encounter     # covariance on the encounter plane 
    p                 # parameter for the poc constraint, dependent on hbr, pcdes, cov_encounter 
    R_tilde           #transformation matrix to transform from ECI frame to the 2D b-plane 
end

# Set of scaling factors to condition an optimization problem.
struct ScalingUnits
    distance_scale
    time_scale
    velocity_scale
    acceleration_scale
end

# struct for optimization problem variables per solve.
struct ProblemVars
    N_period               # number of knot points over the horizon
    nx                     # state dimension
    nu                     # control dimension
    n_days                 # days until TCA
    reference_trajectory   # typically nx × N_period matrix
    all_A                  # vector of A matrices along reference
    all_B                  # vector of B matrices along reference
    dt                     # timestep (scaled units)
    dv_limit               # Δv limit per dt segment (m/s)
end

#Define the scaling units here 
scaling_units = ScalingUnits(100, 1000, 100/1000, 100/(1000^2))



#Define the scaling untis here
# distance_scale = 100
# time_scale = 1000
# velocity_scale = distance_scale/time_scale 
# acceleration_scale = distance_scale/(time_scale^2)

#scaling_units = ScalingUnits(distance_scale, time_scale, velocity_scale, acceleration_scale)

