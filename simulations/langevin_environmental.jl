using Random
using Statistics
using FFTW
using Printf
using NPZ

##############################################################
# PARAMETERS FROM SAD -> RATES
##############################################################

function sad_to_rates(kG, zetaG, alphaG, sigma_e2)
    sigma_d2 = kG * sigma_e2
    mu = 0.5 * zetaG * sigma_d2
    d = 0.5 * (alphaG - 1) * sigma_e2

    d1 = (sigma_d2 + d) / 2
    b1 = (sigma_d2 - d) / 2

    b0 = max(mu, 0.0)
    d0 = max(-mu, 0.0)

    return (
        b0 = b0,
        b1 = b1,
        d0 = d0,
        d1 = d1,
        sigma_d2 = sigma_d2,
        sigma_e2 = sigma_e2,
    )
end

##############################################################
# RNA PARAMETERS
##############################################################

function set_rna_rates(kG, kT; r = 1000.0)
    gamma = kG * r / kT
    return (r = r, gamma = gamma)
end

##############################################################
# LANGEVIN SIMULATION
##############################################################

function simulate_xy(
    x0,
    y0,
    b0,
    b1,
    d0,
    d1,
    sigma_d2,
    sigma_e2,
    r,
    gamma;
    dt,
    T,
)
    counter = 0

    x = zeros(Float64, T)
    y = zeros(Float64, T)

    x[1] = x0
    y[1] = y0

    sqrt_dt = sqrt(dt)

    @inbounds for t in 2:T
        xtm1 = x[t - 1]
        ytm1 = y[t - 1]

        mu_x = (b1 - d1) * xtm1 + (b0 - d0)
        mu_y = r * xtm1 - gamma * ytm1

        noise_x = sqrt(max(sigma_d2 * xtm1 + sigma_e2 * xtm1^2, 0.0))
        noise_y = sqrt(max(r * xtm1 + gamma * ytm1, 0.0))

        dx = mu_x * dt + noise_x * randn() * sqrt_dt
        dy = mu_y * dt + noise_y * randn() * sqrt_dt

        new_x = xtm1 + dx
        new_y = ytm1 + dy

        if new_x < 0.0
            new_x = 0.0
            counter += 1
        end
        if new_y < 0.0
            new_y = 0.0
            counter += 1
        end

        x[t] = new_x
        y[t] = new_y
    end

    return x, y, counter
end

##############################################################
# AUTOCORRELATION
##############################################################

function fftconvolve_full(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    n = length(a) + length(b) - 1
    nfft = nextpow(2, n)

    apad = zeros(Float64, nfft)
    bpad = zeros(Float64, nfft)
    @inbounds apad[1:length(a)] .= a
    @inbounds bpad[1:length(b)] .= b

    fa = rfft(apad)
    fb = rfft(bpad)
    conv = irfft(fa .* fb, nfft)

    return conv[1:n]
end

function autocorr_fft(x::AbstractVector{<:Real})
    xc = Float64.(x) .- mean(x)

    if maximum(abs.(xc)) <= 1e-12
        return ones(Float64, 1)
    end

    result = fftconvolve_full(xc, reverse(xc))
    result = result[(length(result) ÷ 2 + 1):end]

    if result[1] == 0.0
        return ones(Float64, 1)
    end

    return result ./ result[1]
end

function decorrelation_time(x::AbstractVector{<:Real}; threshold = 0.01)
    ac = autocorr_fft(x)

    for (idx, v) in pairs(ac)
        if v < threshold
            return idx - 1
        end
    end

    return length(ac)
end

function decorrelation_time_chunked(
    x::AbstractVector{<:Real};
    threshold = 0.01,
    chunk_size = Int(2e6),
    min_chunk_size = 4096,
    aggregate_quantile = 0.75,
    max_chunks = 20,
    saturation_threshold = 0.25,
    max_chunk_growth_steps = 3,
    max_chunk_size = nothing,
    return_diagnostics = false,
)
    n = length(x)
    chunk_size = max(Int(chunk_size), 1)
    max_chunks = max(Int(max_chunks), 1)
    min_chunk_size = max(Int(min_chunk_size), 1)
    max_chunk_growth_steps = max(Int(max_chunk_growth_steps), 0)

    if max_chunk_size === nothing
        max_chunk_size = n
    else
        max_chunk_size = max(Int(max_chunk_size), 1)
    end

    max_chunk_size = min(max_chunk_size, n)
    effective_chunk_size = min(chunk_size, max_chunk_size)

    if n <= effective_chunk_size
        tau = max(decorrelation_time(x; threshold = threshold), 1)
        if !return_diagnostics
            return tau
        end
        diagnostics = Dict{String, Any}(
            "chunk_size_initial" => chunk_size,
            "chunk_size_used" => n,
            "growth_steps" => 0,
            "n_chunks_used" => 1,
            "saturation_fraction" => Float64(tau >= n),
            "adaptation_exhausted" => false,
        )
        return tau, diagnostics
    end

    growth_steps = 0
    tau_estimate = 1
    saturation_fraction = 0.0
    n_chunks_used = 0
    adaptation_exhausted = false

    while true
        starts = collect(1:effective_chunk_size:n)
        if length(starts) > max_chunks
            select_idx = floor.(Int, range(0, stop = length(starts) - 1, length = max_chunks)) .+ 1
            starts = starts[select_idx]
        end

        taus = Int[]
        saturated_chunks = 0

        for start_idx in starts
            stop_idx = min(start_idx + effective_chunk_size - 1, n)
            chunk = @view x[start_idx:stop_idx]
            if length(chunk) < min_chunk_size
                continue
            end

            tau_chunk = max(decorrelation_time(chunk; threshold = threshold), 1)
            push!(taus, tau_chunk)

            if tau_chunk >= length(chunk)
                saturated_chunks += 1
            end
        end

        n_chunks_used = length(taus)
        if n_chunks_used == 0
            tau_estimate = 1
            saturation_fraction = 0.0
        else
            tau_estimate = max(Int(ceil(quantile(taus, aggregate_quantile))), 1)
            saturation_fraction = saturated_chunks / n_chunks_used
        end

        can_grow = (
            growth_steps < max_chunk_growth_steps &&
            effective_chunk_size < max_chunk_size &&
            effective_chunk_size < n
        )

        if saturation_fraction > saturation_threshold && can_grow
            next_chunk_size = min(effective_chunk_size * 2, max_chunk_size, n)
            if next_chunk_size == effective_chunk_size
                break
            end
            effective_chunk_size = next_chunk_size
            growth_steps += 1
            continue
        end

        adaptation_exhausted = (saturation_fraction > saturation_threshold) && (!can_grow)
        break
    end

    if !return_diagnostics
        return tau_estimate
    end

    diagnostics = Dict{String, Any}(
        "chunk_size_initial" => chunk_size,
        "chunk_size_used" => effective_chunk_size,
        "growth_steps" => growth_steps,
        "n_chunks_used" => n_chunks_used,
        "saturation_fraction" => saturation_fraction,
        "adaptation_exhausted" => adaptation_exhausted,
    )

    return tau_estimate, diagnostics
end

##############################################################
# OUTPUT
##############################################################

function format_scientific_array(values; precision = 3)
    fmt = Printf.Format("%." * string(precision) * "e")
    return map(v -> parse(Float64, Printf.format(fmt, Float64(v))), values)
end

function build_output_filename(params::Dict{String, Any}; chain = 0, outdir = nothing)
    if outdir === nothing
        outdir = joinpath(@__DIR__, "data", "langevin_environmental")
    end

    mkpath(outdir)

    ordered_keys = [
        "kG",
        "zetaG",
        "alphaG",
        "kT",
        "sigma_e2",
        "b1_target",
        "r",
        "x0",
        "y0",
        "T",
        "dt",
        "decor_thresh",
    ]

    pieces = String[]
    for key in ordered_keys
        value = params[key]
        if key == "T"
            push!(pieces, "$(key)=$(Int(value))")
        elseif key == "decor_thresh"
            push!(pieces, @sprintf("%s=%.2e", key, Float64(value)))
        else
            push!(pieces, @sprintf("%s=%.4f", key, Float64(value)))
        end
    end

    push!(pieces, "chain=$(Int(chain))")
    return joinpath(outdir, join(pieces, "__") * ".npz")
end

##############################################################
# MAIN
##############################################################

function parse_main_args(args)
    n_chains = 1
    base_seed = nothing
    T_override = nothing
    dt_override = nothing
    burnin_frac_override = nothing
    burnin_tau_mult_override = nothing
    burnin_min_points_override = nothing

    for arg in args
        if startswith(arg, "--chains=") || startswith(arg, "--n-chains=")
            n_chains = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--seed=")
            base_seed = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--T=")
            T_raw = parse(Float64, split(arg, "=", limit = 2)[2])
            if !isfinite(T_raw) || T_raw < 1
                error("T must be a finite number >= 1")
            end
            T_override = max(Int(round(T_raw)), 1)
        elseif startswith(arg, "--dt=")
            dt_raw = parse(Float64, split(arg, "=", limit = 2)[2])
            if !isfinite(dt_raw) || dt_raw <= 0
                error("dt must be a finite number > 0")
            end
            dt_override = dt_raw
        elseif startswith(arg, "--burnin-frac=")
            burnin_frac_raw = parse(Float64, split(arg, "=", limit = 2)[2])
            if !isfinite(burnin_frac_raw) || burnin_frac_raw < 0 || burnin_frac_raw >= 1
                error("burnin-frac must be in [0, 1)")
            end
            burnin_frac_override = burnin_frac_raw
        elseif startswith(arg, "--burnin-tau-mult=")
            burnin_tau_mult_raw = parse(Int, split(arg, "=", limit = 2)[2])
            if burnin_tau_mult_raw < 0
                error("burnin-tau-mult must be >= 0")
            end
            burnin_tau_mult_override = burnin_tau_mult_raw
        elseif startswith(arg, "--burnin-min-points=")
            burnin_min_points_raw = parse(Int, split(arg, "=", limit = 2)[2])
            if burnin_min_points_raw < 1
                error("burnin-min-points must be >= 1")
            end
            burnin_min_points_override = burnin_min_points_raw
        elseif occursin(r"^[0-9]+$", arg)
            n_chains = parse(Int, arg)
        else
            error(
                "Unknown argument: $(arg). Use --chains=N [--seed=S] [--T=steps] [--dt=value] " *
                "[--burnin-frac=F] [--burnin-tau-mult=M] [--burnin-min-points=P].",
            )
        end
    end

    if n_chains < 1
        error("n_chains must be >= 1")
    end

    return (
        n_chains,
        base_seed,
        T_override,
        dt_override,
        burnin_frac_override,
        burnin_tau_mult_override,
        burnin_min_points_override,
    )
end

function main()
    ##########################################################
    # TARGET PARAMETERS
    ##########################################################

    x0 = 1.0
    y0 = 1.0
    kG = 20.0
    zetaG = 2.0
    alphaG = 1.1
    kT = 50.0
    sigma_e2 = 1.0
    b1_target = 1.0
    r = 5000.0
    dt = 1e-5
    T = Int(5e8)
    decor_thresh = 0.01
    decor_chunk_size = Int(2e6)
    decor_aggregate_quantile = 0.75
    decor_max_chunks = 20
    decor_saturation_threshold = 0.25
    decor_max_chunk_growth_steps = 3
    burnin_frac = 0.1
    burnin_tau_mult = 20
    burnin_min_points = 5000

    n_chains, base_seed, T_override, dt_override, burnin_frac_override, burnin_tau_mult_override, burnin_min_points_override =
        parse_main_args(ARGS)
    if T_override !== nothing
        T = T_override
    end
    if dt_override !== nothing
        dt = dt_override
    end
    if burnin_frac_override !== nothing
        burnin_frac = burnin_frac_override
    end
    if burnin_tau_mult_override !== nothing
        burnin_tau_mult = burnin_tau_mult_override
    end
    if burnin_min_points_override !== nothing
        burnin_min_points = burnin_min_points_override
    end

    ##########################################################
    # BUILD RATES
    ##########################################################

    rates = sad_to_rates(kG, zetaG, alphaG, sigma_e2)
    rna = set_rna_rates(kG, kT; r = r)

    println("Rates:")
    println(rates)
    println(rna)
    println("Number of chains: $(n_chains)")
    println("Simulation settings: T=$(T), dt=$(dt)")
    println(
        "Burn-in settings: burnin_frac=$(burnin_frac), " *
        "burnin_tau_mult=$(burnin_tau_mult), burnin_min_points=$(burnin_min_points)",
    )
    if base_seed !== nothing
        println("Base RNG seed: $(base_seed)")
    end

    params = Dict{String, Any}(
        "x0" => x0,
        "y0" => y0,
        "kG" => kG,
        "zetaG" => zetaG,
        "alphaG" => alphaG,
        "kT" => kT,
        "sigma_e2" => sigma_e2,
        "b1_target" => b1_target,
        "r" => r,
        "dt" => dt,
        "T" => T,
        "decor_thresh" => decor_thresh,
        "decor_chunk_size" => decor_chunk_size,
        "decor_aggregate_quantile" => decor_aggregate_quantile,
        "decor_max_chunks" => decor_max_chunks,
        "decor_saturation_threshold" => decor_saturation_threshold,
        "decor_max_chunk_growth_steps" => decor_max_chunk_growth_steps,
        "burnin_frac" => burnin_frac,
        "burnin_tau_mult" => burnin_tau_mult,
        "burnin_min_points" => burnin_min_points,
    )

    for chain in 0:(n_chains - 1)
        println("=== Chain $(chain + 1)/$(n_chains) (id=$(chain)) ===")
        if base_seed !== nothing
            Random.seed!(base_seed + chain)
        end

        ##########################################################
        # SIMULATION
        ##########################################################

        x, y, counts_boundary = simulate_xy(
            x0,
            y0,
            rates.b0,
            rates.b1,
            rates.d0,
            rates.d1,
            rates.sigma_d2,
            rates.sigma_e2,
            rna.r,
            rna.gamma;
            dt = dt,
            T = T,
        )
        println("Simulation completed. Boundary hits: $(counts_boundary)")

        ##########################################################
        # DECORRELATION
        ##########################################################

        x_tau_steps_pre, x_decor_diag_pre = decorrelation_time_chunked(
            x;
            threshold = decor_thresh,
            chunk_size = decor_chunk_size,
            aggregate_quantile = decor_aggregate_quantile,
            max_chunks = decor_max_chunks,
            saturation_threshold = decor_saturation_threshold,
            max_chunk_growth_steps = decor_max_chunk_growth_steps,
            return_diagnostics = true,
        )
        y_tau_steps_pre, y_decor_diag_pre = decorrelation_time_chunked(
            y;
            threshold = decor_thresh,
            chunk_size = decor_chunk_size,
            aggregate_quantile = decor_aggregate_quantile,
            max_chunks = decor_max_chunks,
            saturation_threshold = decor_saturation_threshold,
            max_chunk_growth_steps = decor_max_chunk_growth_steps,
            return_diagnostics = true,
        )
        max_tau_steps_pre = max(max(x_tau_steps_pre, y_tau_steps_pre), 1)

        n_points = length(x)
        burn_n_from_frac = Int(round(burnin_frac * n_points))
        burn_n_from_tau = burnin_tau_mult * max_tau_steps_pre
        burn_n_discarded = max(burn_n_from_frac, burn_n_from_tau)
        max_discard_allowed = max(n_points - burnin_min_points, 0)
        burn_n_discarded = min(burn_n_discarded, max_discard_allowed)
        burn_start_idx = burn_n_discarded + 1

        x_post = x[burn_start_idx:end]
        y_post = y[burn_start_idx:end]

        x_tau_steps, x_decor_diag = decorrelation_time_chunked(
            x_post;
            threshold = decor_thresh,
            chunk_size = decor_chunk_size,
            aggregate_quantile = decor_aggregate_quantile,
            max_chunks = decor_max_chunks,
            saturation_threshold = decor_saturation_threshold,
            max_chunk_growth_steps = decor_max_chunk_growth_steps,
            return_diagnostics = true,
        )
        y_tau_steps, y_decor_diag = decorrelation_time_chunked(
            y_post;
            threshold = decor_thresh,
            chunk_size = decor_chunk_size,
            aggregate_quantile = decor_aggregate_quantile,
            max_chunks = decor_max_chunks,
            saturation_threshold = decor_saturation_threshold,
            max_chunk_growth_steps = decor_max_chunk_growth_steps,
            return_diagnostics = true,
        )
        max_tau_steps = max(max(x_tau_steps, y_tau_steps), 1)

        x_uncorrelated = x_post[1:max_tau_steps:end]
        y_uncorrelated = y_post[1:max_tau_steps:end]

        println("decorrelation step: $(max_tau_steps)")
        println("decorrelation step pre-burn: $(max_tau_steps_pre)")
        println("burn-in discarded points: $(burn_n_discarded) / $(n_points)")
        println("decorrelation chunk size: $(decor_chunk_size)")
        println("decorrelation aggregation quantile: $(decor_aggregate_quantile)")
        println("decorrelation max chunks: $(decor_max_chunks)")
        println("decorrelation saturation threshold: $(decor_saturation_threshold)")
        println("decorrelation max chunk growth steps: $(decor_max_chunk_growth_steps)")
        x_chunk_size_used = x_decor_diag["chunk_size_used"]
        x_growth_steps = x_decor_diag["growth_steps"]
        x_sat_fraction = Float64(x_decor_diag["saturation_fraction"])
        x_n_chunks_used = x_decor_diag["n_chunks_used"]
        y_chunk_size_used = y_decor_diag["chunk_size_used"]
        y_growth_steps = y_decor_diag["growth_steps"]
        y_sat_fraction = Float64(y_decor_diag["saturation_fraction"])
        y_n_chunks_used = y_decor_diag["n_chunks_used"]
        println(
            "x decor diagnostics: " *
            "chunk_size_used=$(x_chunk_size_used), " *
            "growth_steps=$(x_growth_steps), " *
            @sprintf("saturation_fraction=%.2f, ", x_sat_fraction) *
            "n_chunks_used=$(x_n_chunks_used)",
        )
        println(
            "y decor diagnostics: " *
            "chunk_size_used=$(y_chunk_size_used), " *
            "growth_steps=$(y_growth_steps), " *
            @sprintf("saturation_fraction=%.2f, ", y_sat_fraction) *
            "n_chunks_used=$(y_n_chunks_used)",
        )

        if x_decor_diag["adaptation_exhausted"]
            println(
                "WARNING: x decorrelation estimate may be censored." *
                " Increase decor_chunk_size or decor_max_chunk_growth_steps.",
            )
        end

        if y_decor_diag["adaptation_exhausted"]
            println(
                "WARNING: y decorrelation estimate may be censored." *
                " Increase decor_chunk_size or decor_max_chunk_growth_steps.",
            )
        end

        filename = build_output_filename(params; chain = chain)
        x_uncorrelated_sc = format_scientific_array(x_uncorrelated)
        y_uncorrelated_sc = format_scientific_array(y_uncorrelated)

        out = Dict{String, Any}(
            "x_uncorrelated" => x_uncorrelated_sc,
            "y_uncorrelated" => y_uncorrelated_sc,
            "max_tau_steps" => max_tau_steps,
            "max_tau_steps_preburn" => max_tau_steps_pre,
            "x_tau_steps" => x_tau_steps,
            "y_tau_steps" => y_tau_steps,
            "x_tau_steps_preburn" => x_tau_steps_pre,
            "y_tau_steps_preburn" => y_tau_steps_pre,
            "x_decor_chunk_size_used" => x_decor_diag["chunk_size_used"],
            "y_decor_chunk_size_used" => y_decor_diag["chunk_size_used"],
            "x_decor_growth_steps" => x_decor_diag["growth_steps"],
            "y_decor_growth_steps" => y_decor_diag["growth_steps"],
            "x_decor_saturation_fraction" => x_decor_diag["saturation_fraction"],
            "y_decor_saturation_fraction" => y_decor_diag["saturation_fraction"],
            "x_decor_n_chunks_used" => x_decor_diag["n_chunks_used"],
            "y_decor_n_chunks_used" => y_decor_diag["n_chunks_used"],
            "x_decor_adaptation_exhausted" => x_decor_diag["adaptation_exhausted"],
            "y_decor_adaptation_exhausted" => y_decor_diag["adaptation_exhausted"],
            "x_decor_chunk_size_used_preburn" => x_decor_diag_pre["chunk_size_used"],
            "y_decor_chunk_size_used_preburn" => y_decor_diag_pre["chunk_size_used"],
            "x_decor_growth_steps_preburn" => x_decor_diag_pre["growth_steps"],
            "y_decor_growth_steps_preburn" => y_decor_diag_pre["growth_steps"],
            "x_decor_saturation_fraction_preburn" => x_decor_diag_pre["saturation_fraction"],
            "y_decor_saturation_fraction_preburn" => y_decor_diag_pre["saturation_fraction"],
            "x_decor_n_chunks_used_preburn" => x_decor_diag_pre["n_chunks_used"],
            "y_decor_n_chunks_used_preburn" => y_decor_diag_pre["n_chunks_used"],
            "x_decor_adaptation_exhausted_preburn" => x_decor_diag_pre["adaptation_exhausted"],
            "y_decor_adaptation_exhausted_preburn" => y_decor_diag_pre["adaptation_exhausted"],
            "burn_start_idx" => burn_start_idx,
            "burn_n_discarded" => burn_n_discarded,
            "n_points_total" => n_points,
            "n_points_postburn" => length(x_post),
            "chain" => chain,
            "gamma" => rna.gamma,
            "b0" => rates.b0,
            "b1" => rates.b1,
            "d0" => rates.d0,
            "d1" => rates.d1,
            "sigma_d2" => rates.sigma_d2,
            "sigma_e2_eff" => rates.sigma_e2,
        )

        for (k, v) in params
            out[k] = v
        end

        npzwrite(filename, out)
        println("Saved uncorrelated data: $(filename)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
