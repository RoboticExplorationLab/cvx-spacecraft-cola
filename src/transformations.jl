#transformations 

function rotate_cov_rtn_to_eci(cov_rtn, sat_state)

    r = sat_state[1:3]
    v = sat_state[4:6]

    R_hat = r/norm(r)
    N_hat = cross(r, v)/norm(cross(r,v))
    T_hat = cross(N_hat, R_hat)

    R_mat = [R_hat T_hat N_hat]

    #build 6x6 block matrix (diagonal)
    # [R_mat
    #       R_mat]
    M = kron(Matrix(1.0*I,2,2), R_mat)

    cov_eci = M* cov_rtn *M'

    return cov_eci

end


#add in a rotation to the LVLH frame for the controls 
