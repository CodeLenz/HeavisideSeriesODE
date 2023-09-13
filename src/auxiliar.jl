
########################################################################
#                          Initial conditions                          #
########################################################################

# Defines a function to evaluate the integration constants
function CI_HS1(C_bar::AbstractArray{TF}, F11::AbstractArray{TC}, 
                CbF::AbstractArray{TC},n_excitedDOF::Int64, U0::Vector, 
                V0::Vector) where {TF,TC}

    # Calculates C2
    C2 = (C_bar.-(2*F11))\(V0.+(CbF*U0))

    # Calculates C1
    C1 = @. U0-C2

    return C1, C2

end

########################################################################
#                         Coefficients                                 #
########################################################################
function Coefficients_a_HS1(i::Int64,j::Int64,n_times::Int64,dt::TF,
                            times::Ts,loads::Matrix{TF},int_loads::Matrix{TF}) where {Ts,TF}

    
    # We do not want to compute it if i==n_times
    i==n_times && return zero(TF), zero(TF)

    # Load at i+1
    gL = loads[i+1,j]

    # Load at i
    gl = loads[i,j]

    # Time i
    tl = times[i]

    # Time i+1
    tL = times[i+1]

    # Coeficient a1
    a1 = (gL-gl)/dt

    # Integral of g in the interval. 
    integral = int_loads[i,j]

    # Coeficient a0
    a0 = integral/dt - a1*(tL^2 - tl^2)/(2*dt)

    return a0, a1

end
