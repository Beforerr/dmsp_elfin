using NaNMath
import NaNMath as nm
using Statistics: mean
export init_guess

# https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/plot_exp_cutoff_powerlaw.html

"""
f(E) = A * E^(-γ) * exp(-E/E_c)
"""
struct PowerLawExp{T}
    A::T
    γ::T
    E_c::T
end

(f::PowerLawExp)(E) = f.A * E^(-f.γ) * exp(-E / f.E_c)

"""
    fit(PowerLawExp, E, y)

Fit the model f(E) = A * E^(-γ) * exp(-E/E_c) to data (E, y)
by minimizing ∑[ln(yᵢ) - ln f(Eᵢ)]².

Returns a NamedTuple with fields
- `A`: amplitude
- `γ`: power-law index
- `E_c`: cutoff energy
"""
function fit(::Type{<:PowerLawExp}, E, y)
    N = length(E)
    # Log-transform the data and build the design matrix
    yln = log.(y)
    X = hcat(ones(N), -log.(E), -E)

    # Solve the normal equations via backslash \(x, y)
    # For rectangular A the result is the minimum-norm least squares solution computed by a pivoted QR factorization of A and a rank estimate of A based on the R factor.
    θ = X \ yln
    α, γ, φ = θ

    # Back-transform to original parameters
    A = exp(α)
    E_c = 1 / φ
    return PowerLawExp(A, γ, E_c)
end

# https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/plot_smooth_broken_powerlaw.html
struct SmoothBrokenPowerlaw{T}
    A::T
    γ1::T
    γ2::T
    Eb::T
    m::T
end

(f::SmoothBrokenPowerlaw)(E) = sbpl(E, f.A, f.γ1, f.γ2, f.Eb, f.m)

function fit(::Type{<:SmoothBrokenPowerlaw}, E, y; kw...)
    alg = NonlinearSolve.TrustRegion()
    f = log_sbpl_model_sciml
    p = init_guess(f, E, y)
    prob = NonlinearCurveFitProblem(f, p, E, log.(y))
    sol = solve(prob; alg, kw...)
    return SmoothBrokenPowerlaw(sol.u..., 1.0)
end

function log_sbpl(E, A, γ1, γ2, Eb, m)
    x = E / Eb
    return nm.log(A) - γ1 * nm.log(x) + ((γ1 - γ2) / m) * nm.log(1 + x^m)
end

sbpl(args...) = exp(log_sbpl(args...))

function log_sbpl_model(E, p)
    A, γ1, γ2, Eb = p
    m = 1
    return log_sbpl.(E, A, γ1, γ2, Eb, m)
end

include("model.jl")

for model in (:log_plec_model, :log_sbpl_model, :log_two_pop_model, :log_plec_sbpl_model)
    sciml_model = Symbol(model, :_sciml)
    @eval $sciml_model(p, E) = $model(E, p)
    @eval export $model
    @eval export $sciml_model
    @eval init_guess(::typeof($sciml_model), energies, flux) = init_guess($model, energies, flux)
end


using LeastSquaresOptim, NonlinearSolve
using CurveFit

export remove_nan, fit_flux_two_step

# TODO: ignore the last channel of dmsp_flux

function remove_nan(Xs...)
    idx = mapreduce(.&, Xs) do x
        .!isnan.(x)
    end
    return map(x -> x[idx], Xs)
end


function fit_flux_two_step(flux, energies, Emin; kw...)
    # first fit energy below Emin with `PowerLawExp`
    flux_1 = flux[energies .< Emin]
    Es1 = energies[energies .< Emin]
    f1 = fit(PowerLawExp, Es1, flux_1)

    # second fit the remaining flux of energy above Emin with `SmoothBrokenPowerlaw`
    δEs = energies[energies .>= Emin]
    δflux = flux[energies .>= Emin] - f1.(δEs)
    f2 = fit(SmoothBrokenPowerlaw, δEs, δflux)

    flux_1_modeled = f1.(Es1)
    flux_2_modeled = f2.(δEs) .+ f1.(δEs)

    return (Es1, f1), (δEs, f2), vcat(flux_1_modeled, flux_2_modeled)
end

"""
    msd_log(a, b)

Return the mean squared deviation between the log-transformed two arrays: `mean(abs2, log(a) - log(b))`.
"""
function msd_log(a, b)
    return mean(abs2, log.(a) .- log.(b))
end


"""
    fit_flux_two_step(flux, energies; kw...)

Fit two-step model with optimized transition energy Emin.
Uses actual energy channel values as candidates for Emin.

Returns:
- Best fit results with optimal Emin
"""
function fit_flux_two_step(flux, energies; kw...)
    # Use actual energy values as candidates (instrument channels)
    # Filter to reasonable range and ensure enough points on both sides
    sorted_energies = sort(unique(energies))
    Emins = filter(sorted_energies) do Emin
        n_low = sum(energies .< Emin)
        n_high = sum(energies .>= Emin)
        n_low >= 3 && n_high >= 3
    end

    isempty(Emins) && return nothing, nothing, nothing

    best_fit = nothing
    best_score = Inf
    best_Emin = nothing

    for Emin in Emins
        try
            (Es1, f1), (δEs, f2), flux_modeled = fit_flux_two_step(flux, energies, Emin; kw...)
            score = msd_log(flux, flux_modeled)
            if score < best_score
                best_score = score
                best_Emin = Emin
                best_fit = ((Es1, f1), (δEs, f2), flux_modeled)
            end
        catch
            continue
        end
    end
    return best_fit, best_Emin, best_score
end

export fit_row_parameters

energies(x) = parent(x.dims[1].val)

function fit_row_parameters(flux1, flux2; mlat = nothing)
    ff = vcat(parent(flux1), parent(flux2))
    ee = vcat(energies(flux1), energies(flux2))
    ff, ee = remove_nan(ff, ee)
    n_points = length(ff)

    best_fit, Emin, score = fit_flux_two_step(ff, ee)
    
    if isnothing(best_fit)
        success = false
        params = nothing
        flux_modeled = nothing
    else
        (Es1, f1), (δEs, f2), flux_modeled = best_fit
        params = (f1, f2)
        success = true
    end
    return (; success, flux_modeled, params, n_points, Emin, score)
end
