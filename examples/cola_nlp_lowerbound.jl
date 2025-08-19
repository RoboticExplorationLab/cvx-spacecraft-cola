#NLP formulation for Example 2 (with a lower bound)
#does not solve within the maximum iteration count, which means a more informed initial guess is needed
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


#indexing helper function for setting up the IPOPT problem
function create_idx(nx,nu,N)
    # This function creates some useful indexing tools for Z 
    # x_i = Z[idx.x[i]]
    # u_i = Z[idx.u[i]]

    # our Z vector is [x0, u0, x1, u1, …, xN]
    nz = (N-1) * nu + N * nx # length of Z 
    x = [(i - 1) * (nx + nu) .+ (1 : nx) for i = 1:N]
    u = [(i - 1) * (nx + nu) .+ ((nx + 1):(nx + nu)) for i = 1:(N - 1)]
    
    #constraint indexing for the (N-1) dynamics constraints when stacked up
    c = [(i - 1) * (nx) .+ (1 : nx) for i = 1:(N - 1)]
    nc = (N - 1) * nx # (N-1)*nx 

    
    return (nx=nx,nu=nu,N=N,nz=nz, x= x,u = u, c = c, nc=nc)
end


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

#final time
epcf = epc0 + T

#time series
horizon = LinRange(0, T/scaling_units.time_scale, N_period)

#control upper limit in scaled units
#control_upper_limit = (0.04/(horizon[2]*scaling_units.time_scale))/scaling_units.acceleration_scale

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
Pc_des = 8e-6

#build the conjunction struct 
conjunction = construct_conjunction(reference_trajectory, scaling_units, Pc_des)

#build the indices
idx = create_idx(nx, nu, N_period)

#quadratic cost on the controls
function cost(params, Z)

    N = params.N
    J = 0 

    for i=1:(N-1)

        ui = Z[idx.u[i]]

        J += (norm(ui))^2

    end

    return J

end

#linear dynamcis to have a fair comparison with the sdp method
function dynamics_constraints(params, Z)

    N = params.N 

    c_d = zeros(eltype(Z), idx.nc)

    for i=1:N-1

        c_d[idx.c[i]] = all_A[:,:,i]*Z[idx.x[i]] + all_B[:,:,i]*Z[idx.u[i]] - Z[idx.x[i+1]]

    end

    return c_d

end

  
#dynamics constraints here are linearized. 
#however, they can also be the true nonlinear discreteized dynamics
function equality_constraints(params, Z)

    return [
    
    Z[idx.x[1]] - params.xic; 
    dynamics_constraints(params, Z)
    ]

end


function inequality_constraints(params, Z)

    #inequality_c = zeros(eltype(Z), params.N)

    inequality_c = zeros(eltype(Z), 2*params.N-1)

    #one terminal constraint on the state (meet desired pc)
    #N-1 control constraints on an upper bound on the L2 control

    rel_encounter = conjunction.R_tilde*(conjunction.sat1_tca[1:3] + Z[idx.x[end]][1:3] - conjunction.sat2_tca[1:3])
    
    inequality_c[1]= rel_encounter'*inv(conjunction.cov_encounter)*rel_encounter - conjunction.p;

    #upper bound constraint
    for i=1:params.N-1
        inequality_c[i+1] = norm(Z[idx.u[i]]) -  ((0.01/(horizon[2]*scaling_units.time_scale))/scaling_units.acceleration_scale)
    end


    #lower bound constraint
    for i=1:params.N-1

        #max iterations with the squared norm (that is the formulation in the sdp)
        inequality_c[params.N + i] = (norm(Z[idx.u[i]]))^2 - (((0.008/(horizon[2]*scaling_units.time_scale))/scaling_units.acceleration_scale)/5)^2

    end

    return inequality_c

end


#set the parameters 
#xic is the initial condition
params = (N = N_period, xic = [(0.1*ones(3))/scaling_units.distance_scale; (0.01*ones(3))/scaling_units.velocity_scale], idx = idx)

#create initial guess 
function create_initial_guess(idx)

    initial_guess = zeros(idx.nz)

    for i=1:size(idx.x)[1]

        initial_guess[idx.x[i]] = reference_trajectory[:,i]

    end
    #leave the controls initial guess as zeros (not good for new version of ForwardDiff, set small delta)

    # for i=1:size(idx.u)[1]
    #     initial_guess[idx.u[i]] = 1e-10*ones(idx.nu)
    # end

    return initial_guess
end

#initial guess 
x0 = create_initial_guess(idx)

#create the bounds for the equality constraints 
x_l = -Inf * ones(idx.nz)
x_u = Inf * ones(idx.nz)

#create the bounds for the inequality constraints 
c_l = -Inf*ones(2*N_period-1)
c_u = Inf*ones(2*N_period-1)

#update lower bound for the first inequaltiy 
c_l[1] = 0.0

#update the upper bound for the controls constraint 

for i=2:N_period 

    c_u[i] = 0.0

end

#update the bound for the lower bound constraint 

for i=1:N_period-1

    #lower bound is zero, upper bound is good at INf
    c_l[N_period+i] = 0

end

diff_type = :auto 

#solve the problem 

cost(params, x0) 

x = fmincon(cost::Function,
                                  equality_constraints::Function,
                                  inequality_constraints::Function,
                                  x_l::Vector,
                                  x_u::Vector,
                                  c_l::Vector,
                                  c_u::Vector,
                                  x0::Vector,
                                  params::NamedTuple,
                                  diff_type::Symbol;
                                  tol = 1e-7,
                                  c_tol = 1e-7,
                                  max_iters = 1_000,
                                  verbose = true)

#extract the solution
#unpack the solution 
x_sol = zeros(nx, N_period)
u_sol = zeros(nu, N_period-1)


for i = 1:N_period

    x_sol[:,i] = x[idx.x[i]]
end


for i=1:N_period-1

    u_sol[:,i] = x[idx.u[i]]

end

#cost for the solution 

cost_sol = 0

for i=1:N_period-1

    cost_sol += norm(u_sol[:,i])^2
end

cost_sol 

#get the pc ellipse in scaled units 
ellipse_pc_scaled = plot_pc_ellipse(conjunction.p, conjunction.cov_encounter)

#in km
ellipse_pc_SI = ellipse_pc_scaled*scaling_units.distance_scale/1000

#relative position in the b-plane plot
relative_position_encounter = conjunction.R_tilde*(conjunction.sat1_tca[1:3] + x_sol[1:3,end] - conjunction.sat2_tca[1:3])

#get in km 
relative_position_SI = relative_position_encounter*scaling_units.distance_scale/1000

plot(ellipse_pc_SI[1,:], ellipse_pc_SI[2,:], xlabel ="Bx (km)", ylabel="Bz (km)", title="Deviation in Encounter Plane", label=false)
scatter!([relative_position_SI[1]], [relative_position_SI[2]], label="nlp solve")
