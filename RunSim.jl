module RunSim

export run_sim

using DSP


"""
    run_sim(p::Int, T::Int, init::Symbol, update::Symbol, coupling::Bool)

Runs simulation on system size `2^p` for `T` steps. Returns power spectra at
wavenumbers `2^q` for `q=1:p`.

`init` describes the initialization method, which can be `:rand_init` for random
configurations or `:dyck_init` for configurations drawn uniformly from the Dyck paths.

`update` describes the update method, which can be `:open_update` or
`:periodic_update` for open or periodic boundary conditions respectively.

`coupling` determines whether or not updates on parallel simulations are in the
same (randomly chosen) direction, causing simulations to couple
"""
function run_sim(p::Int, T::Int, init::Symbol, update::Symbol, coupling::Bool)

    N = 2^p  # Size of the system

    path = eval(init)(2^(p-1)) # init functions return system of size 2N

    ft_traces = [Complex{Float64}[] for _ in 1:p]

    for t in 1:T

        path = eval(update)(path, coupling) # Don't like this much!!

        if t % N == 0
            # Collect spectra only after all sites updated on average
            ft = fft([count_ones(site) for site in path])
            for q in 1:p
                push!(ft_traces[q],ft[2^q])
            end
        end

    end

    [Periodograms.periodogram(ft_traces[q]).power for q in 1:p]

end


"""
    periodic_update(config, coupling::Bool)

Update the configuration of a periodic Fredkin chain at a randomly chosen bond.

If `coupling=true` the updates are chosen (randomly) to be in the same direction.
"""
function periodic_update(config, coupling::Bool)

    if coupling
        gate = fredkin_gate_coupling
    else
        gate = fredkin_gate
    end


    size = length(config)

    # We do two updates, with the control bit on the left...

    control_left = rand(collect(1:size))

    if control_left == size
        config[1], config[2] = gate(config[size], config[1], config[2])
    elseif control_left == size - 1
        config[size], config[1] = gate(config[size - 1], config[size], config[1])
    else
        config[control_left + 1], config[control_left + 2] =
        gate(config[control_left], config[control_left + 1], config[control_left + 2])
    end

    # ... and on the right

    control_right = rand(collect(1:size))

    if control_right == 1
        config[size - 1], config[size] = gate(config[1], config[size - 1], config[size])
    elseif control_right == 2
        config[size], config[1] = gate(config[2], config[size], config[1])
    else
        config[control_right - 2], config[control_right - 1] =
        gate(config[control_right], config[control_right - 2], config[control_right - 1])
    end


    return config
end


"""
    open_update(config, coupling::Bool)

Update the configuration of a open Fredkin chain at a randomly chosen bond.

If `coupling=true` the updates are chosen (randomly) to be in the same direction.
"""
function open_update(config, coupling)

    if coupling
        gate = fredkin_gate_coupling
    else
        gate = fredkin_gate
    end


    size = length(config)

    # We do two updates, with the control bit on the left...

    control_left = rand(collect(1:size))

    if control_left < size-1
        config[control_left + 1], config[control_left + 2] =
        gate(config[control_left], config[control_left + 1], config[control_left + 2])
    end

    # ... and on the right

    control_right = rand(collect(1:size))

    if control_right > 2
        config[control_right - 2], config[control_right - 1] =
        gate(config[control_right], config[control_right - 2], config[control_right - 1])
    end


    return config
end

"""
    fredkin_gate(C::UInt, I1::UInt, I2::UInt)

Computes the two outputs of a Fredkin gate with inputs `C` (control),
`I1` and `I2`. If `C = 1`, `I1` and `I2` are swapped if different,
otherwise they remain unchanged.
"""
function fredkin_gate(C::UInt, I1::UInt, I2::UInt)
    S = xor(I1, I2) & C
    return xor(I1, S), xor(I2, S)
end

"""
    fredkin_gate_coupling(C::UInt, I1::UInt, I2::UInt)

Computes the two outputs of a Fredkin gate with inputs `C` (control), `I1` and `I2`.
If `C = 1`, `I1` and `I2` are swapped if different in a randomly chosen direction.
Otherwise they remain unchanged.
"""
function fredkin_gate_coupling(C::UInt, I1::UInt, I2::UInt)
    S = rand([I1 & ~I2, ~I1 & I2]) & C
    return xor(I1, S), xor(I2, S)
end

"""
    rand_init(n::Int)

Return 2n random `UInt64`s for random initial conditions
"""
function rand_init(n::Int)
    return rand(UInt64, 2n)
end

"""
    dyck(n::Int)

Returns a Dyck path of 1s and 0s of length 2n sampled uniformly.

From http://bit.ly/2FPeeOd
"""
function dyck(n::Int)
    # Note n is the half the length to ensure Dyck constraint
    path = append!(ones(Int, n),zeros(Int, n+1))
    path = shuffle(path) # Random permutation
    nadir = indmin(cumsum(path))
    circshift(path, -nadir)[1:2n]
end

"""
    init_dyck_paths(n::Int)

Return 2n `UInt64`s that describe 64 Dyck paths sampled uniformly.
"""
function dyck_init(n::Int)

    configs = hcat([dyck(n) for _ in 1:64]...)

    integer_configs = []

    for j in 1:2n
        bit_string = join(configs[j,:])
        append!(integer_configs, parse(UInt64, bit_string, 2)) # Note it's an unsigned integer!
    end

    return integer_configs
end

end
