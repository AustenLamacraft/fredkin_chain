module FredkinChain

export find_power_spectra, lorentzian_time

#include("RunSim.jl")
using RunSim
using CurveFit

"""
    find_power_spectra(N_sim::Int, size_exp::Int, time_steps::Int;
        init::Symbol=:dyck_init, update::Symbol=:open_update, coupling::Bool=false, ncores::Int=8)

Finds the average power spectra over `N_sim` simulations run on `ncores` cores,
for a chain of length `2^(size_exp)`.

`init` describes the initialization method, which can be `:rand_init` for random
configurations or `:dyck_init` for configurations drawn uniformly from the Dyck paths.

`update` describes the update method, which can be `:open_update` or `:periodic_update`
for open or periodic boundary conditions respectively.

`coupling` determines whether or not updates on parallel simulations are in the
same (randomly chosen) direction.
"""
function find_power_spectra(N_sim::Int, size_exp::Int, time_steps::Int; init::Symbol=:dyck_init, update::Symbol=:open_update, coupling::Bool=false, ncores::Int=8)

    init in [:rand_init, :dyck_init] ||
                error("""
                      Unknown initialization: $init
                             Allowed methods:
                               :rand_init
                               :dyck_init""")

    update in [:open_update, :periodic_update] ||
                error("""
                      Unknown update: $update
                             Allowed methods:
                               :open_update
                               :periodic_update""")


    sum_over_cores = @parallel (+) for i=1:ncores
        mean_power_specs = [zeros(time_steps/2^size_exp) for _ in 1:size_exp]
        for i in 1:ceil(Int, N_sim / ncores)
            print("running simulation: ", i, "\n")
            mean_power_specs += run_sim(size_exp, time_steps, init, update, coupling)
        end
        mean_power_specs
    end

    sum_over_cores = sum_over_cores / N_sim
end

"""
    lorentzian_time(x, y)

Fit to Lorentzian plus a constant. Interpret coefficient of `x^2` in denominator as squared time
"""
function lorentzian_time(x, y)
    return sqrt(abs(curve_fit(RationalPoly, x.^2, y, 1, 1).den.a[2]))
end

end
