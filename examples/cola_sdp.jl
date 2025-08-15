#cola sdp 
#activate local environment
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
using MathOptInterface 
using DelimitedFiles

include("../src/dynamics_integrate.jl")
include("../src/conjunction.jl")
include("../src/problem_objects.jl")
include("../src/transformations.jl")

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

#build the conjunction struct 
conjunction = construct_conjunction(reference_trajectory, scaling_units)

#x and u optimization variables
xtraj = Variable(nx,N_period)
utraj = Variable(nu, N_period-1)

#each column is a matrix here from the moment matrix
xxtraj = Variable(nx^2, N_period)
uutraj = Variable(nu^2, N_period-1)
xutraj = Variable(nx*nu, N_period-1)

#psd constraints 
#/succeq is for the psd constraint 
constraints = Constraint[ [1 xtraj[:,1]' utraj[:,1]'; 
                           xtraj[:,1] reshape(xxtraj[:,1], nx, nx) reshape(xutraj[:,1], nx, nu); 
                           utraj[:,1] reshape(xutraj[:,1], nx, nu)' reshape(uutraj[:,1], nu, nu)] ⪰ 0 
                        ];

for k=2:(N_period-1)
    push!(constraints, [1 xtraj[:,k]' utraj[:,k]'; 
                           xtraj[:,k] reshape(xxtraj[:,k], nx, nx) reshape(xutraj[:,k], nx, nu); 
                           utraj[:,k] reshape(xutraj[:,k], nx, nu)' reshape(uutraj[:,k], nu, nu)] ⪰ 0)

end

#terminal psd constraint 
push!(constraints, [1 xtraj[:,N_period]'; 
                    xtraj[:,N_period] reshape(xxtraj[:,N_period], nx, nx)] ⪰ 0)

#dynamics constraints 
for k=1:(N_period-1)

    #discrete dynamics
    push!(constraints, xtraj[:,k+1] == all_A[:,:,k]*xtraj[:,k] + all_B[:,:,k]*utraj[:,k])
    
    #these discrete dynamics on the quadratic terms (to promote tightness)
    push!(constraints, reshape(xxtraj[:,k+1], nx, nx) == all_A[:,:,k]*reshape(xxtraj[:,k], nx, nx)*all_A[:,:,k]' + all_A[:,:,k]*reshape(xutraj[:,k], nx, nu)*all_B[:,:,k]' + (all_A[:,:,k]*reshape(xutraj[:,k], nx, nu)*all_B[:,:,k]')' + all_B[:,:,k]*reshape(uutraj[:,k], nu, nu)*all_B[:,:,k]')

end

#convex poc constraint
push!(constraints, tr(inv(conjunction.cov_encounter)*conjunction.R_tilde*(conjunction.sat1_tca[1:3]*conjunction.sat1_tca[1:3]' + conjunction.sat1_tca[1:3]*xtraj[1:3,N_period]' - conjunction.sat1_tca[1:3]*conjunction.sat2_tca[1:3]' + xtraj[1:3,N_period]*conjunction.sat1_tca[1:3]' + reshape(xxtraj[:,N_period], nx, nx)[1:3, 1:3] - xtraj[1:3,N_period]*conjunction.sat2_tca[1:3]' - conjunction.sat2_tca[1:3]*conjunction.sat1_tca[1:3]' - conjunction.sat2_tca[1:3]*xtraj[1:3, N_period]' + conjunction.sat2_tca[1:3]*conjunction.sat2_tca[1:3]')*conjunction.R_tilde') >= conjunction.p)

for k=1:N_period-1 

    #upper bound for ctrl
    ub = (dv_limit/(horizon[2]*scaling_units.time_scale))/scaling_units.acceleration_scale

    #tight with this upper bound constraint
    push!(constraints, norm(utraj[:,k]) <= ub)

end

#set this as the initial condition away from the reference 
x0 = [(0.1*ones(3))/scaling_units.distance_scale; (0.01*ones(3))/scaling_units.velocity_scale]

push!(constraints, xtraj[:,1] == x0)
push!(constraints, reshape(xxtraj[:,1], nx, nx) == x0*x0')

#minimize the squared terms of the control
problem = minimize(sum(uutraj[1,:]) + sum(uutraj[5,:]) + sum(uutraj[9,:]), constraints); 

optimizer = Mosek.Optimizer()

Convex.solve!(problem, ()-> optimizer, silent=false)

#extract the solution 
xtraj_value = xtraj.value
utraj_value = utraj.value

#each column is a matrix here from the moment matrix. 
xxtraj_value = xxtraj.value
uutraj_value = uutraj.value
xutraj_value = xutraj.value

#subtraction check on the consistency of x_traj and x_traj x_traj'
check_x = zeros(nx,nx, N_period)

for i=1:N_period

    #check that the outer product matrices are consistent 
    check_x[:,:,i] = xtraj_value[:,i]*xtraj_value[:,i]' - reshape(xxtraj_value[:,i], nx, nx)

end 

check_u = zeros(nu,nu, N_period-1)

for i=1:N_period-1

    #check that the outer product matrices are consistent 
    check_u[:,:,i] = utraj_value[:,i]*utraj_value[:,i]' - reshape(uutraj_value[:,i], nu, nu)

end

check_xu = zeros(nx, nu, N_period-1)

for i=1:N_period-1

    check_xu[:,:,i] = xtraj_value[:,i]*utraj_value[:,i]' - reshape(xutraj_value[:,i], nx, nu) 

end

#all close to zero so we're good
sum(check_x)
sum(check_u)
sum(check_xu)

#Now check moment matrix rank to ensure tightness of the problem

all_e = zeros(N_period)

for k=1:N_period-1 

    M = [1 xtraj_value[:,k]' utraj_value[:,k]'; 
                           xtraj_value[:,k] reshape(xxtraj_value[:,k], nx, nx) reshape(xutraj_value[:,k], nx, nu); 
                           utraj_value[:,k] reshape(xutraj_value[:,k], nx, nu)' reshape(uutraj_value[:,k], nu, nu)] 

    eig_vals = eigen(M).values 

    all_e[k] = eig_vals[end]/eig_vals[end-1]


end

M_N = [1 xtraj_value[:,N_period]'; 
                    xtraj_value[:,N_period] reshape(xxtraj_value[:,N_period], nx, nx)]

eig_vals_N = eigen(M_N).values 
all_e[N_period] = eig_vals_N[end]/eig_vals_N[end-1]


plot(all_e, xlabel="Knot Point", ylabel="Eigenvalue Ratio", title="Tightness Check", label="eigen value ratio",  yscale=:log10, ylim= (10e4, 10e8)) 


#Plot the results

#plot the deviation from the reference trajectory
plot(horizon[1:N_period]*scaling_units.time_scale/T, xtraj_value[1,:]*scaling_units.distance_scale/1000, title="Deviation from Reference", xlabel = "Orbit Revs", ylabel="Distance (km)", label="dx")
plot!(horizon[1:N_period]*scaling_units.time_scale/T, xtraj_value[2,:]*scaling_units.distance_scale/1000, label="dy")
plot!(horizon[1:N_period]*scaling_units.time_scale/T, xtraj_value[3,:]*scaling_units.distance_scale/1000, label="dz")

#plot the controls in mm/s2
plot(horizon[1:N_period-1]*scaling_units.time_scale/T, utraj_value[1,:]*scaling_units.acceleration_scale*1000, title="SDP Method Maneuvers", label = "ax")
plot!(horizon[1:N_period-1]*scaling_units.time_scale/T, utraj_value[2,:]*scaling_units.acceleration_scale*1000, label="ay")
plot!(horizon[1:N_period-1]*scaling_units.time_scale/T, utraj_value[3,:]*scaling_units.acceleration_scale*1000, label="az", xlabel = "Orbit Revs", ylabel="acceleration (mm/s2)")


#plot the maneuver magnitude
normz = zeros(N_period-1)

for i=1:N_period-1

    normz[i] = norm(utraj_value[:,i])

end
plot(horizon[1:N_period-1]*scaling_units.time_scale/T, normz*scaling_units.acceleration_scale*1000, title="Maneuver Magnitude", xlabel = "Orbit Revs", ylabel = "Magnitude (mm/s2)", legend=false)

#get the pc ellipse in scaled units 
ellipse_pc_scaled = plot_pc_ellipse(conjunction.p, conjunction.cov_encounter)


#plot the manevuer effect in the encounter plane (b-plane)
relative_position_encounter = conjunction.R_tilde*(conjunction.sat1_tca[1:3] + xtraj_value[1:3,N_period] - conjunction.sat2_tca[1:3])
relative_position_encounter_no_maneuver = conjunction.R_tilde*(conjunction.sat1_tca[1:3] - conjunction.sat2_tca[1:3])

#where is this delta on the encounter plane...
plot(ellipse_pc_scaled[1,:]*scaling_units.distance_scale/1000, ellipse_pc_scaled[2,:]*scaling_units.distance_scale/1000, title="SDP Method Maneuver Effect", xlabel = "Ex (km)", ylabel="Ez (km)", label="pc ellipse 1e-6")
scatter!([relative_position_encounter[1]*scaling_units.distance_scale/1000], [relative_position_encounter[2]*scaling_units.distance_scale/1000], label="with maneuver")
scatter!([relative_position_encounter_no_maneuver[1]*scaling_units.distance_scale/1000], [relative_position_encounter_no_maneuver[2]*scaling_units.distance_scale/1000], label="no maneuver")
scatter!([0], [0], label="spacecraft 2")

#how good is the linearization 
#here we run the controls on the true nonlinear dynamics model
sim_reference = zeros(9, N_period)
sim_reference[:,1] = [(reference_trajectory[:,1] + xtraj_value[:,1]); utraj_value[:,1]] 
for i=1:N_period-2

    #simulate the reference with added control 
    solution_one_dt = just_dynamics_integrate(sim_reference[:,i], horizon[i], horizon[i+1])

    #get the state at the next timestep 
    next_state = get_state(solution_one_dt, N_period)[:,end]

    sim_reference[:, i+1] = next_state

    #replace the control with the control from the opt ctrl problem 
    sim_reference[7:9, i+1] = utraj_value[:, i+1]

end

#integrate the final state
solution_one_dt = just_dynamics_integrate(sim_reference[:,N_period-1], horizon[N_period-1], horizon[N_period])

final_state = get_state(solution_one_dt, N_period)[:,end]

sim_reference[:, N_period] = final_state

true_deltas = sim_reference[1:6, :] - reference_trajectory
linearized_deltas = xtraj_value

#what is the difference between the true deltas and the linearized deltas
linearization_error = true_deltas - linearized_deltas 

#cm level eror on positon, and mm/s level error on velocity. linearization is pretty good
linearization_error_SI_position = linearization_error[1:3, :]*scaling_units.distance_scale
linearization_error_SI_velocity = linearization_error[4:6, :]*scaling_units.velocity_scale