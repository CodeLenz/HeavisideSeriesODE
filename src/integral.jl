#
# Evaluates the integral of the loads using the trapezoidal
# method
#
function Integral_load(loads::Matrix{TF},dt::TF) where TF

    # Initializes the matrix
    int_loads = zero(loads)

    # Loop over times
    for i = 1:size(loads,1)-1
        int_loads[i,:] .= 0.5*dt*(loads[i+1,:].+loads[i,:])
    end #i

    # Return the integrals
    return int_loads

end
