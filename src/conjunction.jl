include("transformations.jl")
include("problem_objects.jl")

#define the hard body radius (m). same for both examples
hbr = 10

#projection onto the b-plane matrix
E = [1.0 0.0 0.0; 0.0 0.0 1.0]

#define conjunction and the bplane coordinate system
#get the miss vector from an arbitrary cdm and apply it to our system
object_1 = 
        [
            -2269.916517,
            -6492.472918,
            66.637818,
            -1.376274071,
            0.567532572,
            7.552308423
        ]*1e3


object_2 = 
        [
            -2269.004683,
            -6492.879724,
            66.63485, 
            -1.427895217,
            0.411252835,
            -7.551930038 
        ]*1e3
        
#get this delta to make up a state for satellite 2
delta = object_2 - object_1

tca_epoch = Epoch("2024-06-17T17:41:37.496000Z")

#CDM Covariances 

#object 1 position uncertainty in RTN frame 
CR_R                               =4655.970725060175        #[m**2]
CT_R                               =-505376.5903145122       #[m**2]
CT_T                               =57273213.98347315        #[m**2]
CN_R                               =156.5352780669326        #[m**2]
CN_T                               =-16184.77656518719       #[m**2]
CN_N                               =105.6790981099849        #[m**2]

#object 1 velocity uncertainty in RTN frame 
CRDOT_R                            =558.4656768050252        #[m**2/s]
CRDOT_T                            =-63279.59560893906       #[m**2/s]
CRDOT_N                            =17.88525706193112        #[m**2/s]
CRDOT_RDOT                         =69.91596689472058        #[m**2/s**2]
CTDOT_R                            =-1.847136248141935       #[m**2/s]
CTDOT_T                            =190.2775260014937        #[m**2/s]
CTDOT_N                            =-0.07860641018958604     #[m**2/s]
CTDOT_RDOT                         =-0.2103419520082558      #[m**2/s**2]
CTDOT_TDOT                         =0.0009038026907550972    #[m**2/s**2]
CNDOT_R                            =0.03277294155665673      #[m**2/s]
CNDOT_T                            =-2.350226224721231       #[m**2/s]
CNDOT_N                            =0.05078815288763895      #[m**2/s]
CNDOT_RDOT                         =0.002598039703224794     #[m**2/s**2]
CNDOT_TDOT                         =-0.00004842123406964361  #[m**2/s**2]
CNDOT_NDOT                         =0.0001028765153830954    #[m**2/s**2]

#object 2 position uncertainty in RTN frame 
CR_R2                               =37.35063992668302        #[m**2]
CT_R2                               =-8960.511954729523       #[m**2]
CT_T2                               =3039967.700836536        #[m**2]
CN_R2                               =-2.076044837796598       #[m**2]
CN_T2                               =688.9207383826407        #[m**2]
CN_N2                               =16.503780525201          #[m**2]

#object 2 velocity uncertainty in RTN frame 
CRDOT_R2                            =9.949906086495822        #[m**2/s]
CRDOT_T2                            =-3373.748961190803       #[m**2/s]
CRDOT_N2                            =-0.7645809588466025      #[m**2/s]
CRDOT_RDOT2                         =3.744182297782079        #[m**2/s**2]
CTDOT_R2                           =-0.01250291409863463     #[m**2/s]
CTDOT_T2                            =3.462959330715508        #[m**2/s]
CTDOT_N2                            =0.0007910083097934923    #[m**2/s]
CTDOT_RDOT2                         =-0.003844387843034295    #[m**2/s**2]
CTDOT_TDOT2                         =0.000004923687437770908  #[m**2/s**2]
CNDOT_R2                            =-0.00130214592391184     #[m**2/s]
CNDOT_T2                            =0.3510639256634225       #[m**2/s]
CNDOT_N2                            =-0.002902558806392447    #[m**2/s]
CNDOT_RDOT2                         =-0.0003896927195015376   #[m**2/s**2]
CNDOT_TDOT2                         =0.0000004417957899003916 #[m**2/s**2]
CNDOT_NDOT2                         =0.00001891244436111274   #[m**2/s**2]

object_1_cov_rtn = [
    CR_R     CT_R     CN_R     CRDOT_R     CTDOT_R     CNDOT_R;
    CT_R     CT_T     CN_T     CRDOT_T     CTDOT_T     CNDOT_T;
    CN_R     CN_T     CN_N     CRDOT_N     CTDOT_N     CNDOT_N;
    CRDOT_R  CRDOT_T  CRDOT_N  CRDOT_RDOT  CTDOT_RDOT  CNDOT_RDOT;
    CTDOT_R  CTDOT_T  CTDOT_N  CTDOT_RDOT  CTDOT_TDOT  CNDOT_TDOT;
    CNDOT_R  CNDOT_T  CNDOT_N  CNDOT_RDOT  CNDOT_TDOT  CNDOT_NDOT
]

# covariance in the RTN frame (already in meters and seconds)
object_2_cov_rtn = [
    CR_R2     CT_R2     CN_R2     CRDOT_R2     CTDOT_R2     CNDOT_R2;
    CT_R2     CT_T2     CN_T2     CRDOT_T2     CTDOT_T2     CNDOT_T2;
    CN_R2     CN_T2     CN_N2     CRDOT_N2     CTDOT_N2     CNDOT_N2;
    CRDOT_R2  CRDOT_T2  CRDOT_N2  CRDOT_RDOT2  CTDOT_RDOT2  CNDOT_RDOT2;
    CTDOT_R2  CTDOT_T2  CTDOT_N2  CTDOT_RDOT2  CTDOT_TDOT2  CNDOT_TDOT2;
    CNDOT_R2  CNDOT_T2  CNDOT_N2  CNDOT_RDOT2  CNDOT_TDOT2  CNDOT_NDOT2
]

function construct_conjunction(reference_trajectory, scaling_units, Pc_des)

    sat1_tca = reference_trajectory[:,end]

    sat2_tca = sat1_tca + [delta[1:3]./scaling_units.distance_scale; delta[4:6]./scaling_units.velocity_scale]

    object_1_cov_eci = rotate_cov_rtn_to_eci(object_1_cov_rtn, sat1_tca)  
    object_2_cov_eci = rotate_cov_rtn_to_eci(object_2_cov_rtn, sat2_tca)

    r_relative = sat1_tca[1:3,end] - sat2_tca[1:3]
    v_relative = sat1_tca[4:6,end] - sat2_tca[4:6]


    #define the encounter plane axes to rotate eci into this
    y_axis = v_relative/norm(v_relative)

    z_axis = cross(r_relative, v_relative)/norm(cross(r_relative, v_relative))

    x_axis = cross(y_axis, z_axis)

    eci2encounter = [x_axis'; y_axis'; z_axis']

    #combined covariance in the eci frame
    cov_combined_eci = object_1_cov_eci[1:3, 1:3] + object_2_cov_eci[1:3, 1:3]

    cov_combined_encounter = eci2encounter* cov_combined_eci* eci2encounter'
    
    #get a 2x2 covariance matrix in the encounter plane 
    #assumption that there is no uncertainty in the relative direction. 3D problem becomes 2D 

    #equivalent notation in the paper
    cov_encounter = E*cov_combined_encounter*E'

    #replaced c_ellipsoid_test with C

    p = log(hbr^4/ (Pc_des^2*4*det(cov_encounter)))

    #the units are meters squared in cov_encoutner, scale to custom units
    cov_encounter_scaled = cov_encounter./(scaling_units.distance_scale)^2

    #convert to encounter plane and project onto b-plane
    R_tilde = E*eci2encounter

    #create the conjunction object and return it 
    conjunction = ConjunctionEvent(sat1_tca, sat2_tca, cov_encounter_scaled, p, R_tilde)

    return conjunction 

end

#this pc ellipse will be in scaled units
function plot_pc_ellipse(c_ellipse, C_ca)

    #number of points for plotting ellipse 
    N_ellipse = 100

    t = range(0, 2*pi, length=N_ellipse)

    #eigen decomposition of the covariance matrix 
    eig = eigen(C_ca)
    λ = eig.values 
    V = eig.vectors 

    #build the scaled ellipse 
    circle = [cos.(t)'; sin.(t)']
    A =  sqrt(c_ellipse)*Diagonal(sqrt.(λ))
    ellipse = V*A*circle #rotation then scale 

    #x_ellipse = ellipse

    return ellipse

end