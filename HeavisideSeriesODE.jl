__precompile__()
module HeavisideSeriesODE
    
    # If possible, set optimization to 3
    if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
        @eval Base.Experimental.@optlevel 3
    end

    # load packages
    using LinearAlgebra
    
    # Load the local subroutines
    include("auxiliar.jl")
    include("integral.jl")
    include("conjugated.jl")
    include("notconjugated.jl")
    include("solver.jl")

    # Export the methods
    export Solve_HS1, Integral_load


end
