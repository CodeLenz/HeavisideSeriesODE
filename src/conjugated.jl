# Defines a function to evaluate the HS response when CbF and F11 are 
# complex conjugate and the initial conditions are not homogeneous
function CbFFconj!(response::Array{TC}, K::AbstractArray, M::AbstractArray, 
                   CbF::AbstractArray, F11::AbstractArray, F11_squared::AbstractArray,
                   expF11_delta::AbstractArray, dim::Integer,
                   n_excitedDOF::Integer, times::Ts, n_times::Integer, input_vectors::Matrix{TF}, 
                   loads::Matrix{TF},int_loads::Matrix{TF},C2::Vector) where {Ts,TC,TF}
 
     # Calculates matrix Omega_2
     Ω_2 = K*CbF
 
     # Calculates matrix Omega_4
     Ω_4 = K.-(M*F11_squared)
 
     # Calculates matrix Omega_5
     Ω_5 = Ω_4*F11
 
     # Initializes and calculates a matrix for the vector x_ij for each
     # matrix Gamma_i 
     X_1 = Array{TC}(undef, dim, n_excitedDOF)
     X_2 = Array{TC}(undef, dim, n_excitedDOF)
     X_4 = Array{TC}(undef, dim, n_excitedDOF)
     X_5 = Array{TC}(undef, dim, n_excitedDOF)
 
     # Initializes the vector of gamma vectors, which convey information 
     # from one iteration to the next
 
     # gamma1
     Γ_1 = zeros(TC, dim, n_excitedDOF)
 
     # Initializes the vectors of sums ch
     χ_0 = zeros(TC, n_excitedDOF)
     χ_1 = zeros(TC, n_excitedDOF)
 
     # Factorizes omegas # TODO -> use LinearSolve
     lu_Ω_1 = cholesky(Symmetric(K))
     lu_Ω_2 = lu!(Ω_2)
     lu_Ω_4 = lu!(Ω_4)
     lu_Ω_5 = lu!(Ω_5)
 
     # Time step
     Δt = times[2]-times[1]
 
     # Iterates through the degrees of freedom
     for j=1:n_excitedDOF
 
         # Gamma_1
         X_1[:,j] .= lu_Ω_1\input_vectors[:,j]
 
         # Gamma_2
         X_2[:,j] .= lu_Ω_2\input_vectors[:,j]
 
         # Gamma_4
         X_4[:,j] .= lu_Ω_4\input_vectors[:,j]
 
         # Gamma_5
         X_5[:,j] .= lu_Ω_5\input_vectors[:,j]
 
         # Calculates the coefficients of HS
         a0, a1 = Coefficients_a_HS1(1,j,n_times,Δt,times,loads,int_loads)
 
         # Adds the vector of initial conditions directly at gamma1 to 
         # spare computational cost
         Γ_1[:,j] = @. ((a1*X_5[:,j])-(((times[1]*a1)+a0)*X_4[:,j])+C2)
 
         # Sum chi_0
         χ_0[j] += a0
 
         # Sum chi_1
         χ_1[j] += a1
 
         # Adds the first response point
         @. response[:,1] = (response[:,1]+(((times[1]*X_1[:,j])-
                             (2*real.(X_2[:,j])))*χ_1[j])+(X_1[:,j]*χ_0[j])+
                             (2*(real.(Γ_1[:,j]))))
 
     end #j
 
     # Pre-alocate gamma_mult1j
     Γ_mult1j = Vector{TC}(undef,dim)

     # Iterates through the time points
     for i=2:n_times
 
        # Current time point
        t_current = times[i]
 
        # Iterates through the degrees of freedom
        for j=1:n_excitedDOF
 
            # Calculates the coefficients of HS
            a0, a1 = Coefficients_a_HS1(i,j,n_times,Δt,times,loads,int_loads)
    
            # Calculates the coefficients
            cji0 = a0-χ_0[j]
            cji1 = a1-χ_1[j]
    
            # Multiplies the gamma vectors by the exponentials of the delta
            Γ_mult1j .= expF11_delta*Γ_1[:,j]
    
            # Calculates the response and adds the homogeneous response
            @. response[:,i] = (response[:,i]+(2*real.(Γ_mult1j))+(((t_current*X_1[:,j])-
                                (2*real.(X_2[:,j])))*χ_1[j])+(χ_0[j]*X_1[:,j]))
    
            # Updates the gamma vector
            Γ_1[:,j] = @. (Γ_mult1j+(cji1*X_5[:,j])-(((cji1*t_current)+cji0)*X_4[:,j]))
    
            # Updates the chi summations
            χ_0[j] += cji0
            χ_1[j] += cji1
    
         end  # j
  
     end # i
    
 end
 
