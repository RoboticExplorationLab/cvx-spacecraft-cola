#cola lineariztion approach 
using Pkg 
Pkg.activate(joinpath(@__DIR__, ".." ))

using LinearAlgebra
using SatelliteDynamics
using Plots 
using DifferentialEquations
using ForwardDiff
using Convex
using Mosek
using MosekTools
using DelimitedFiles

include("../src/dynamics_integrate.jl")
include("../src/conjunction.jl")
include("../src/problem_objects.jl")
include("../src/transformations.jl")
include("../src/fmincon.jl")


# initial time for sim
epc0 = Epoch(2012, 11, 8, 12, 0, 0, 0.0)

#Orbit we want (around ISS orbit conditions)
iss1 = [6871e3, 0.00064, 51.6450, 1.0804, 27.7899, 190]; 

# Convert osculating elements to Cartesean state
# returns position and velocity (m, m/s). This is the intial position
eci0_1 = sOSCtoCART(iss1, use_degrees=true)

#number of revs to avoid the conjunction
n_revs = 1

T = orbit_period(iss1[1])*n_revs

#50 knot points per number of revs
N_period = 50*n_revs

#state size
nx = 6 

#control size 
nu = 3

#set the delta v limit per timestep 
dv_limit = 0.04 #m/s

#final time
epcf = epc0 + T

#time series
horizon = LinRange(0, T/scaling_units.time_scale, N_period)

#initial state
eci0_1_scaled = [eci0_1[1:3]./scaling_units.distance_scale; eci0_1[4:6]./scaling_units.velocity_scale; zeros(3)]

#integrate a whole period  
solution_just_dynamics = just_dynamics_integrate(eci0_1_scaled, 0.0, T/scaling_units.time_scale)

#extract the state
all_states_just_dynamics  = get_state(solution_just_dynamics, N_period) 

reference_trajectory = zeros(nx, N_period)

controls_trajectory = zeros(nu, N_period)

for k=1:N_period

    #save the reference trajectory at the dedicated timesteps
    reference_trajectory[:,k] = solution_just_dynamics(horizon).u[k][1:6]

    controls_trajectory[:,k] = solution_just_dynamics(horizon).u[k][7:9]

end

#get the discrete dynamics jacobian at each timestep
all_A = zeros(nx, nx, N_period-1)
all_B = zeros(nx, nu, N_period-1)
 
for i=1:N_period-1

    all_A[:,:,i] =  ForwardDiff.jacobian(dx->just_dynamics_integrate(dx, horizon[i], horizon[i+1]).u[end][1:6], [reference_trajectory[:,i]; zeros(3)])[1:6, 1:6]
    all_B[:,:,i] =  ForwardDiff.jacobian(dx->just_dynamics_integrate(dx, horizon[i], horizon[i+1]).u[end][1:6], [reference_trajectory[:,i]; zeros(3)])[1:6, 7:9]

end

#set the desired Pc 
Pc_des = 1e-6

#build the conjunction struct 
conjunction = construct_conjunction(reference_trajectory, scaling_units, Pc_des)

#get the pc ellipse in scaled units 
ellipse_pc_scaled = plot_pc_ellipse(conjunction.p, conjunction.cov_encounter)


#here the decision variables are delta v and delta x 
#chain all the delta



#optimization variables
Δx_i = Variable(nx, N_period) 
u_i = Variable(nu, N_period-1)

#paramters based on the linearization along the ellipse
ellipse_a = Variable((2,1))
ellipse_b = Variable(1)

#match the sdp initial condition
cons = Constraint[Δx_i[:,1] == [(0.1*ones(3))/scaling_units.distance_scale; (0.01*ones(3))/scaling_units.velocity_scale]]

#dynamics constraint 
for k=1:N_period-1 

    push!(cons, zeros(6) == all_A[:,:,k]*Δx_i[:,k] + all_B[:,:,k]*u_i[:,k] - Δx_i[:,k+1])

end

#relative position on the b-plane
relative_position_encounter_i = (conjunction.R_tilde*((conjunction.sat1_tca[1:3] + Δx_i[1:3, end]) - conjunction.sat2_tca[1:3]))

#half-place pc approx. constraint
cons_pc = push!(cons, ellipse_a'*relative_position_encounter_i - ellipse_b >= 0)

obj = 0 

#cost function
for k=1:N_period-1

    obj += sumsquares(u_i[:,k])

end

for k=1:N_period-1 

    push!(cons, norm(u_i[:,k]) <=(0.04/(horizon[2]*scaling_units.time_scale))/scaling_units.acceleration_scale)

end


#create the problem object
prob = minimize(obj, cons)


#linearize around the entire ellipse and find the lowest acceleration (delta v / t)
#make a contour plot of lowest delta v
all_a = zeros(2, size(ellipse_pc_scaled)[2])
all_b = zeros(size(ellipse_pc_scaled)[2])
all_maneuvers = zeros(3, N_period-1, size(ellipse_pc_scaled)[2])
all_delta_x = zeros(nx, N_period, size(ellipse_pc_scaled)[2])
all_norms_SI = zeros(size(ellipse_pc_scaled)[2])
all_status = zeros(size(ellipse_pc_scaled)[2])



#solve the problem with all possible linearizations on the ellipsoid 
for i=1:size(ellipse_pc_scaled)[2]

    println("solving iteration: ", i)
    Δr_ca_test_1 = ellipse_pc_scaled[:,i]
    ∇f_1 = 2*inv(conjunction.cov_encounter)* Δr_ca_test_1

    #calculate a and b 
    all_a[:,i] = ∇f_1
    all_b[i] = dot(∇f_1, Δr_ca_test_1)

    set_value!(ellipse_a, all_a[:,i])
    set_value!(ellipse_b, all_b[i])

    fix!(ellipse_a)
    fix!(ellipse_b)

    Convex.solve!(prob, ()-> Mosek.Optimizer(), silent=true)

    free!(ellipse_a)
    free!(ellipse_b)

    if Δx_i.value == nothing || u_i.value == nothing
        all_delta_x[:,:,i] = zeros(nx, N_period)
        all_maneuvers[:,:,i] = zeros(nu, N_period-1)
        all_status[i] = Int(prob.status)
    else
        all_delta_x[:,:,i] = Δx_i.value
        all_maneuvers[:,:,i] = u_i.value
        all_status[i] = Int(prob.status)
    end

end


#for each of the 100 control trajectories, get the norm
all_normz = zeros(size(ellipse_pc_scaled)[2])

#L2 norm or sum? 
for i=1:size(ellipse_pc_scaled)[2]

    maneuver_norm = 0

    for k=1:N_period-1 

        maneuver_norm += sum(abs.(all_maneuvers[:,k,i]))

    end

    all_normz[i] = maneuver_norm 

end


#filter all the nonzero all normz 
#this is because some instances do not solve because the 
#half plane is too conservative 
all_normz = filter(x -> x != 0, all_normz)

minimum_accel = findmin(all_normz)

minimum_index = minimum_accel[2]
best_maneuvers_SI = all_maneuvers[:,:, minimum_index]*scaling_units.acceleration_scale
best_delta_x = all_delta_x[:,:, minimum_index]

plot(best_maneuvers_SI[1,:])
plot!(best_maneuvers_SI[2,:])
plot!(best_maneuvers_SI[3,:])


#in KM 
plot(ellipse_pc_scaled[1,:]*scaling_units.distance_scale/1000, ellipse_pc_scaled[2,:]*scaling_units.distance_scale/1000, label="Pc contour 1e-6")

#Plot the tangent line: a₁ x + a₂ y = b → y = (b - a₁ x) / a₂
x_vals = range(-50, 100; length=100)
if abs(all_a[2, minimum_index]) > 1e-8
    y_vals = (all_b[minimum_index] .- all_a[1, minimum_index] .* x_vals) ./ all_a[2, minimum_index]
    plot!(x_vals*scaling_units.distance_scale/1000, y_vals*scaling_units.distance_scale/1000, label="Best linearization", color=:red)
else
    # vertical line (a[2] ≈ 0)
    x_line = fill(all_b[:,minimum_index] / all_a[1,minimum_index], 100)
    y_line = range(-4, 4; length=100)
    plot!(x_line*scaling_units.distance_scale/1000, y_line*scaling_units.distance_scale/1000, label="Best Linearization", color=:red)
end

#plot the relative change
final_best_delta = best_delta_x[:,end]

encounter_delta = conjunction.R_tilde*((conjunction.sat1_tca[1:3] + final_best_delta[1:3]) - conjunction.sat2_tca[1:3])

scatter!([encounter_delta[1]]*scaling_units.distance_scale/1000, [encounter_delta[2]]*scaling_units.distance_scale/1000, label="best delta (cvx)")