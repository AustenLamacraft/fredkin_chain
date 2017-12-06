module FredkinChain

export find_power_spectra

# Natural thing is to have a type called FredkinChain, for which we can specify the intial conditions and update.

# Update can be periodic or open
# Init can be random or Dyck

# One option is to pass these as Init or Update types


# Example from https://github.com/JuliaOpt/GLPKMathProgInterface.jl/blob/master/src/GLPKInterfaceLP.jl


# type GLPKSolverLP <: AbstractMathProgSolver
#     presolve::Bool
#     method::Symbol
#     opts
#     function GLPKSolverLP(;presolve::Bool=false, method::Symbol=:Simplex, opts...)
#         method in [:Simplex, :Exact, :InteriorPoint] ||
#             error("""
#                   Unknown method for GLPK LP solver: $method
#                          Allowed methods:
#                            :Simplex
#                            :Exact
#                            :InteriorPoint""")
#         new(presolve, method, opts)
#     end
# end

# Have a chain type with methods:
# Initialize
# Run the system
# Functions running simulation must be available on all workers
# Advantage is we can have a different chain for each kind of simulation
# Different module for running simulation that will be available to all

# Write a function that takes a chain and returns power spectrum

#include("RunSim.jl")
using RunSim

"""
    find_power_spectra(N_sim::Int, size_exp::Int, time_steps::Int; ncores::Int=8)

Finds the average power spectra over `N_sim` simulations for a chain of length `2^(size_exp)`
"""
function find_power_spectra(N_sim::Int, size_exp::Int, time_steps::Int; ncores::Int=8)

    sum_over_cores = @parallel (+) for i=1:ncores
        mean_power_specs = [zeros(time_steps/2^size_exp) for _ in 1:size_exp]
        for i in 1:ceil(Int, N_sim / ncores)
            print("running simulation: ", i, "\n")
            mean_power_specs += run_sim(size_exp, time_steps)
        end
        mean_power_specs
    end

    sum_over_cores = sum_over_cores / N_sim
end

end
