#dynamics and integrate methods
include("problem_objects.jl")

#two body continous dynamics model in scaled units
function two_body_dynamics_scaled(xu_scaled, time_s)

    x_scaled = xu_scaled[1:6]
    u_scaled = xu_scaled[7:9]

    x_original = [x_scaled[1:3].*scaling_units.distance_scale; x_scaled[4:6].*scaling_units.velocity_scale]

    t_original  = time_s*scaling_units.time_scale

    #add the elapsed time to the initial epoch 
    epc = epc0 + t_original 
    
    r = x_original[1:3] #satellite position in inertial frame
    v = x_original[4:6] #satellite velocity in inertial frame
        
    PN = bias_precession_nutation(epc)
    Earth_r = earth_rotation(epc)
    rpm  = polar_motion(epc) 

    R = rpm*Earth_r*PN
    
    #Compute the sun and moon positions in ECI frame
    r_sun = sun_position(epc)
    r_moon = moon_position(epc)
    
    #define the acceleration variable
    a = zeros(eltype(x_scaled), 3)
    
    #compute acceleration caused by Earth gravity
    #modeled by a spherical harmonic gravity field
    n_grav = 10
    m_grav = 10
    #main contribution in acceleration
    a+= accel_gravity(x_original, R, n_grav, m_grav)
    
    #atmospheric drag
    #compute the atmospheric density from density harris priester model
    ρ = density_harris_priester(epc,r)
    
    #computes acceleration due to drag in inertial directions
    cd = 2.0 #drag coefficient
    area_drag = 0.1 #in m2 #area normal to the velocity direction
    m = 1.0
    
    a += accel_drag(x_original, ρ, m, area_drag, cd, Array{Real,2}(PN))
        
    area_srp = 1.0
    coef_srp = 1.8

    a_srp = accel_srp(x_original, r_sun, m, area_srp, coef_srp)

    #acceleration due to external bodies
    a_sun = accel_thirdbody_sun(x_original, r_sun)
    
    a_moon = accel_thirdbody_moon(x_original, r_moon)

    #add in external bodies acceleration and srp 
    a += a_srp + a_sun + a_moon 
          
    #scale back to custom units 
    xdot_scaled = x_original[4:6] ./ scaling_units.velocity_scale

    #assume additive acceleration on the continous dynamics as our control input
    vdot_scaled = (a ./ scaling_units.acceleration_scale) + u_scaled

    #zero order hold on the controls
    udot_scaled = zeros(3)
    
    x_dot_scaled = [xdot_scaled; vdot_scaled; udot_scaled]

    return x_dot_scaled
    
end


#Methods to integrate the dynamics and get the solutions

#x_0[1:6] -> position, velocity 
#x_0[7:9] -> control (acceleration)
function just_dynamics_integrate(x_0, t1, t2)
    
    tspan = (t1, t2)
    prob = ODEProblem(justdynamics_DFJL!, x_0, tspan)
    sol = solve(prob, TsitPap8(), abstol=1e-12, reltol=1e-12)
    
    return sol
    
end

#solve using DifferentialEquations.jl
function justdynamics_DFJL!(du, u, p, t)
    
    du[1:9] = two_body_dynamics_scaled(u[1:9], t)
    
end

#gets the state and for the entire solution
function get_state(solution)
    
    N = size(solution.u)[1]

    all_states = zeros(9, N)

    for i=1:N
        all_states[:,i] = solution.u[i][1:9]
    end
    
    #all states and all stm are functions of t
    #solution.t is the time
    return all_states
end






# function get_jacobians(N_period, nx, nu)

#     all_A = zeros(nx, nx, N_period-1)
#     all_B = zeros(nx, nu, N_period-1)

#     for i=1:N_period-1

#         all_A[:,:,i] =  ForwardDiff.jacobian(dx->just_dynamics_integrate(dx, horizon[i], horizon[i+1]).u[end][1:6], [reference_trajectory[:,i]; zeros(3)])[1:6, 1:6]
#         all_B[:,:,i] =  ForwardDiff.jacobian(dx->just_dynamics_integrate(dx, horizon[i], horizon[i+1]).u[end][1:6], [reference_trajectory[:,i]; zeros(3)])[1:6, 7:9]

#     end
    

# end