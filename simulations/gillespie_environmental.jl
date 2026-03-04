using Random
using Statistics
using FFTW
using Printf
using NPZ

##############################################################
# GILLESPIE SIMULATION
##############################################################

function gillespie_gene_expression(
    t_max,
    x0,
    y0,
    b0,
    b1,
    d1,
    sigma,
    r,
    gamma;
    max_events = Int(1e7),
)
    if any(p -> p < 0, (t_max, b0, b1, d1, sigma, r, gamma))
        error("All rates/parameters must be non-negative.")
    end
    if x0 < 0 || y0 < 0
        error("Initial conditions must be non-negative integers.")
    end
    if max_events < 1
        error("max_events must be >= 1")
    end

    t = 0.0
    x = Int(x0)
    y = Int(y0)

    times = Float64[t]
    xs = Int[x]
    ys = Int[y]

    half_sig2 = 0.5 * sigma * sigma
    event_count = 0
    hit_max_events = false

    while t < t_max
        if event_count >= max_events
            hit_max_events = true
            break
        end

        b_x = b1 + half_sig2 * x
        d_x = d1 + half_sig2 * x

        a1 = b_x * x
        a2 = b0
        a3 = d_x * x
        a4 = r * x
        a5 = gamma * y

        a0 = a1 + a2 + a3 + a4 + a5
        if a0 <= 0.0
            break
        end

        tau = randexp() / a0
        t_next = t + tau
        if t_next > t_max
            break
        end

        u = rand() * a0
        c1 = a1
        c2 = c1 + a2
        c3 = c2 + a3
        c4 = c3 + a4

        if u < c1
            x += 1
        elseif u < c2
            x += 1
        elseif u < c3
            if x > 0
                x -= 1
            end
        elseif u < c4
            y += 1
        else
            if y > 0
                y -= 1
            end
        end

        t = t_next
        event_count += 1
        push!(times, t)
        push!(xs, x)
        push!(ys, y)
    end

    return times, xs, ys, event_count, hit_max_events
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
        outdir = joinpath(@__DIR__, "data", "gillespie_environmental")
    end

    mkpath(outdir)

    ordered_keys = [
        "b0",
        "b1",
        "d1",
        "sigma",
        "r",
        "gamma",
        "x0",
        "y0",
        "t_max",
        "decor_thresh",
    ]

    pieces = String[]
    for key in ordered_keys
        value = params[key]
        if key == "x0" || key == "y0"
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
# CLI
##############################################################

function parse_main_args(args)
    n_chains = 1
    base_seed = nothing
    t_max_override = nothing
    max_events_override = nothing
    burnin_frac_override = nothing
    burnin_tau_mult_override = nothing
    burnin_min_points_override = nothing

    for arg in args
        if startswith(arg, "--chains=") || startswith(arg, "--n-chains=")
            n_chains = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--seed=")
            base_seed = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--tmax=") || startswith(arg, "--t-max=") || startswith(arg, "--T=")
            t_raw = parse(Float64, split(arg, "=", limit = 2)[2])
            if !isfinite(t_raw) || t_raw <= 0
                error("t_max must be a finite number > 0")
            end
            t_max_override = t_raw
        elseif startswith(arg, "--max-events=")
            m_raw = parse(Int, split(arg, "=", limit = 2)[2])
            if m_raw < 1
                error("max-events must be >= 1")
            end
            max_events_override = m_raw
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
                "Unknown argument: $(arg). Use --chains=N [--seed=S] [--tmax=TIME] [--max-events=M] " *
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
        t_max_override,
        max_events_override,
        burnin_frac_override,
        burnin_tau_mult_override,
        burnin_min_points_override,
    )
end

##############################################################
# MAIN
##############################################################

function main()
    # Gillespie baseline params (same logic as the Python prototype)
    t_max = 20000.0
    x0 = 10
    y0 = 0
    b0 = 0.5
    b1 = 0.8
    d1 = 1.0
    sigma = 1.0
    r = 100.0
    gamma = 50.0
    max_events = Int(2e8)

    # Decorrelation settings
    decor_thresh = 0.01
    decor_chunk_size = Int(2e6)
    decor_aggregate_quantile = 0.75
    decor_max_chunks = 20
    decor_saturation_threshold = 0.25
    decor_max_chunk_growth_steps = 3
    burnin_frac = 0.1
    burnin_tau_mult = 20
    burnin_min_points = 5000

    n_chains, base_seed, t_max_override, max_events_override, burnin_frac_override, burnin_tau_mult_override, burnin_min_points_override =
        parse_main_args(ARGS)
    if t_max_override !== nothing
        t_max = t_max_override
    end
    if max_events_override !== nothing
        max_events = max_events_override
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

    println("Parameters:")
    println((
        b0 = b0,
        b1 = b1,
        d1 = d1,
        sigma = sigma,
        r = r,
        gamma = gamma,
        x0 = x0,
        y0 = y0,
        t_max = t_max,
        max_events = max_events,
    ))
    println("Number of chains: $(n_chains)")
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
        "b0" => b0,
        "b1" => b1,
        "d1" => d1,
        "sigma" => sigma,
        "r" => r,
        "gamma" => gamma,
        "t_max" => t_max,
        "max_events" => max_events,
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

        times, x, y, n_events, hit_max_events = gillespie_gene_expression(
            t_max,
            x0,
            y0,
            b0,
            b1,
            d1,
            sigma,
            r,
            gamma;
            max_events = max_events,
        )

        println("Simulation completed. Events: $(n_events). Final t: $(times[end])")
        if hit_max_events
            println("WARNING: reached max_events=$(max_events) before t_max=$(t_max)")
        end

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
        times_post = times[burn_start_idx:end]

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
        times_uncorrelated = times_post[1:max_tau_steps:end]

        println("decorrelation step: $(max_tau_steps)")
        println("decorrelation step pre-burn: $(max_tau_steps_pre)")
        println("burn-in discarded points: $(burn_n_discarded) / $(n_points)")

        filename = build_output_filename(params; chain = chain)

        out = Dict{String, Any}(
            "times_uncorrelated" => times_uncorrelated,
            "x_uncorrelated" => format_scientific_array(x_uncorrelated),
            "y_uncorrelated" => format_scientific_array(y_uncorrelated),
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
            "burn_time" => times[burn_start_idx],
            "n_points_total" => n_points,
            "n_points_postburn" => length(x_post),
            "n_events" => n_events,
            "final_time" => times[end],
            "hit_max_events" => hit_max_events,
            "chain" => chain,
            "b0" => b0,
            "b1" => b1,
            "d1" => d1,
            "sigma" => sigma,
            "r" => r,
            "gamma" => gamma,
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
