# HeavisideSeriesODE

This repository contains the computer implementation of the solution procedures developed in **to be included** 
for solving coupled systems of second order ODEs with constant coefficients 

 $M A(t) + C V(t) + K Y(t) = F(t)$

with initial conditions

 $Y(0) = U0$

and

 $V(0) = V0$

where $t$ is the independent variable, $Y(t)$ a $n \times 1$ vector (dependent variables), $V$ its first derivative with respect to $t$ and $A$ its second derivative. Matrices  $M$, $C$ and $K$ are $n \times n$. Vector $F(t)$ is addressed in the following.
 
As this package is not yet registered, you must install it by using

```julia
using Pkg
Pkg.add("https://github.com/CodeLenz/HeavisideSeriesODE.jl.git")
```

The main method is

```julia
y = Solve_HS1(M::AbstractArray{TF}, C::AbstractArray{TF}, K::AbstractArray{TF}, 
              times::Ts, loads::Matrix{TF}, dofs_loads::Matrix{TF}; 
              U0=[], V0=[], tol_conj=1E-6, int_loads::Union{Matrix{TF}, Nothing}=nothing )
```

where $times_{nt}$ is a StepRange or a vector with the $nt$ discrete times, $loads_{ndof \times nload}$ is a matrix containing the load values for each excited DOF (collumns) for each discrete time (rows). It is assumed that the time step between two consecutive times is always the same.

# Example

Simple $3$ DOFs system with a senoidal excitation. 

```julia
using HeavisideSeriesODE

M = [5.0 0.0 0.0 ;
     0.0 2.0 0.0 ;
     0.0 0.0 3.0]  


K = [6.0 -4.0  0.0 ;
    -4.0 10.0 -3.0 ;
     0.0 -3.0  7.0]*1E2

C =  K*1E-6

# Assuming a time span of 10 seconds
times = range(start=0.0,stop=10.0,length=100)

# The load is applied at the second DOF
dofs_loads = zeros(3,1);
dofs_loads[2,1] = 1.0

# a simple sin with constant frequency
loads = zeros(length(times),1);
loads[:,1] .= sin.(10.0 * times)

# Evaluate the response. Each row is a DOF
# and each column the complete solution at 
# discrete time point
y = Solve_HS1(M,C,K,times,loads,dofs_loads)
```
