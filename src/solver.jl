"""
Evaluates the response of 

MA(t) + CV(t) + KU(t) = F(t)

using the First Order Heaviside Series. 

Inputs:


    M,C,K:: square coefficient matrices of the same size. The size is internaly used
            as dim

    times:: a vector or a steprange with the discrete times. It is assumed that those
            times are equally spaced in time (constant time step). The length of this
            collection is internaly used as n_times.


    loads::  matrix with n_loads columns (one for each loaded DOF) and n_times rows.

    dofs_loads:: matrix with dim rows and n_loads columns. This matrix has 1 in the 
                position of the loaded dof. It is assumed that the sequence of columns
                respect the order assumed in matrix loads.

    Named arguments:

    U0, V0:: vectors with initial conditions. If empty, we assume homogeneous (0) conditions

    tol_conj:: tolerance to check for conjugacy 

    int_loads:: matrix with n_loads columns (one for each loaded DOF) and n_times rows. It
                contains the integral of each load (column) at each interval. If not informed
                the integrals are approximated using the trapezoidal rule

    Output:

    response: complex matrix with U(t) for each DOF and discrete time. The rows correspond to each
            DOF and the collumns to each discrete time

"""
function Solve_HS1(M::AbstractArray{TF}, C::AbstractArray{TF}, K::AbstractArray{TF}, 
                   times::Ts, loads::Matrix{TF}, dofs_loads::Matrix{TF}; 
                   U0=[], V0=[], tol_conj=1E-6, int_loads::Union{Matrix{TF}, Nothing}=nothing ) where {Ts,TF}


    ############################################
    # INITIAL VERIFICATIONS AND INITIALIZATIONS         
    ############################################

    # Number of time points
    n_times = length(times)

    # Assert that we have at least 2 time points
    n_times >=2 || throw("HS1:: we need at least two time points")

    # Number of excited dofs
    n_excitedDOF = size(dofs_loads,2)

    # Assert that we have at least one column (excited dof)
    n_excitedDOF >=1 || throw("HS1:: we need at least one excited DOF")

    # Check if matrices are consistent 
    size(K)==size(M)==size(C) || throw("HS1:: Matrices M,K and C must have the same dimensions")

    # Problem's dimension
    dim,plim = size(M)

    # Check if matrices are square
    dim==plim || throw("HS1:: M,C and K must be square matrices")
    
    # Assert that loads is consistent 
    size(loads)==(n_times,n_excitedDOF) || throw("HS1:: loads size is not consistent") 

    # Assert that dofs_loads is consistent
    size(dofs_loads,1)==dim || throw("HS1:: input_vectors is not consistent with the number of DOFs") 

    # Assumes that the time discretization is constant
    Δt = times[2]-times[1]

    # If int_loads is empty, we can evaluate by using the trapezoidal rule
    if isnothing(int_loads)

        # Trapezoidal rule
        int_load = Integral_load(loads,Δt) 

    else
        
        # Otherwise, check if the size is OK 
        size(int_loads)==(n_times,n_excitedDOF) || throw("HS1:: dimensions of int_load are not consistent") 

        # And use an internal alias
        int_load = int_loads
        
    end # if isempty

    # Check initial condition U0 and evaluated its norm
    if isempty(U0)
        resize!(U0,dim)
        fill!(U0,zero(TF))
        norm_U0 = zero(TF)
    else
        length(U0)==dim || throw("HS1:: size of U0 is not consistent")
        norm_U0 = norm(U0)
    end


    # Check initial condition V0 and evaluated its norm
    if isempty(V0)
        resize!(V0,dim)
        fill!(V0,zero(TF))
        norm_V0 = zero(TF)
    
    else
        length(V0)==dim || throw("HS1:: size of V0 is not consistent")
        norm_V0 = norm(V0)
    end

    ############################################
    # INITIALIZE THE MATRICES        
    ############################################

    # Factorizes the mass matrix
    cholM = cholesky(Symmetric(M))

    # Calculates C_bar and K_bar
    # and convert to arrays as Julia does not evaluate
    # sqrt and exp of sparse arrays
    C_bar = Array(cholM\C)
    K_bar = Array(cholM\K)

    # Square of C_bar. 
    C_bar2 = C_bar*C_bar

    # Calculates F211
    F11 = 0.5*(C_bar .+ sqrt(C_bar2.- 4*K_bar))
    
    # Initializes the response matrix
    response = zeros(ComplexF64, dim, n_times)

    # Calculates F11 squared
    F11_squared = F11*F11

    # Calculates exp(-delta t*F11)
    expF11_delta = exp(-Δt*F11)
    
    # Calculates CbF 
    CbF = C_bar .- F11

    ############################################
    # INTEGRATION CONSTANTS
    ############################################
    if norm(U0)!=0 || norm(V0)!=0

        # Calculates initial conditions
        C1, C2 = CI_HS1(C_bar, F11, CbF, n_excitedDOF, U0, V0)

    else

        # Homogeneous initial conditions
        C1 = zeros(TF,dim)
        C2 = zeros(TF,dim)

    end

    # Compute the norm of the conjugacy
    norm_conjugacy = norm((CbF).-conj(F11))

    ############################################
    # CALL THE APROPRIATE KERNEL         
    ############################################
    if norm_conjugacy<tol_conj

        # Evaluate the response
        CbFFconj!(response,K,M,CbF,F11,F11_squared,expF11_delta,dim,
                n_excitedDOF,times,n_times,dofs_loads,loads,int_load,C2)

    else

        # Calculates the exponential of Cb-F11
        expCF_delta = exp(-Δt*CbF)

        # Evaluate the response        
        CbFFnotconj!(response,K,M,CbF,F11,F11_squared,expF11_delta,expCF_delta,
                    dim,n_excitedDOF,times,n_times,dofs_loads,loads,int_load,C1,C2)

    end


    # Return the response U(t)
    return real.(response)

end
