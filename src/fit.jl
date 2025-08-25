using NaNMath
import NaNMath as nm
using Statistics: mean
export init_guess, TwoStepModel, PowerLaw, SmoothBrokenPowerlaw

# https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/index.html
# https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/plot_exp_cutoff_powerlaw.html
# https://fjebaker.github.io/SpectralFitting.jl/dev/

include("model.jl")

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

"""
    fit(PowerLaw, E, y)

Fit the model f(E) = A * E^(-γ) to data (E, y)
by minimizing ∑[ln(yᵢ) - ln f(Eᵢ)]².

Returns a NamedTuple with fields
- `A`: amplitude
- `γ`: power-law index
"""
function fit(::Type{<:PowerLaw}, E, y)
    N = length(E)
    # Log-transform the data and build the design matrix
    yln = log.(y)
    X = hcat(ones(N), -log.(E))

    # Solve the normal equations
    θ = X \ yln
    α, γ = θ

    # Back-transform to original parameters
    A = exp(α)
    return PowerLaw(A, γ)
end

function fit(::Type{<:SmoothBrokenPowerlaw}, E, y; kw...)
    alg = NonlinearSolve.TrustRegion()
    f = log_sbpl_model_sciml
    p = init_guess(f, E, y)
    prob = NonlinearCurveFitProblem(f, p, E, log.(y))
    sol = solve(prob; alg, kw...)
    return SmoothBrokenPowerlaw(sol.u..., 1.0)
end

function log_sbpl_model(E, p)
    A, γ1, γ2, Eb = p
    m = 1
    return log_sbpl.(E, A, γ1, γ2, Eb, m)
end


function log_plec_model(E, p)
    A, γ, Ec = p
    return nm.log(A) .- γ .* nm.log.(E) .- (E ./ Ec)
end

function log_powerlaw_model(E, p)
    A, γ = p
    return nm.log(A) .- γ .* nm.log.(E)
end

function init_guess(::typeof(log_plec_model), energies, flux)
    i = argmax(flux)
    γ0 = 0.5
    E0 = energies[i]
    Ec0 = energies[end]
    A = flux[i] / E0^γ0 * exp(E0 / Ec0)
    return [A, γ0, Ec0]
end

function init_guess(::typeof(log_sbpl_model), energies, flux)
    i = argmax(flux)
    m = 1
    γ1, γ2 = -5.0, 3.0
    E0 = energies[i]
    Eb = 1.0e3
    A = flux[i] / exp(log_sbpl_model(E0, [1, γ1, γ2, Eb]))
    return [A, γ1, γ2, Eb]
end

for model in (:log_plec_model, :log_sbpl_model)
    sciml_model = Symbol(model, :_sciml)
    @eval $sciml_model(p, E) = $model(E, p)
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


function fit_flux_two_step(flux, energies, Emin; model1 = PowerLawExp, model2 = SmoothBrokenPowerlaw, kw...)
    # first fit energy below Emin with model1 (default: PowerLawExp)
    flux_1 = flux[energies .<= Emin]
    Es1 = energies[energies .<= Emin]
    f1 = fit(model1, Es1, flux_1)

    # second fit the remaining flux of energy above Emin with model2 (default: SmoothBrokenPowerlaw)
    δEs = energies[energies .> Emin]
    δflux = flux[energies .> Emin] - f1.(δEs)
    f2 = fit(model2, δEs, δflux)

    # Create combined model
    model = TwoStepModel(f1, f2, Emin)
    flux_modeled = model.(energies)

    return model, flux_modeled
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
- TwoStepModel: Combined model structure
- best_Emin: Optimal transition energy
- best_score: Fit quality score (lower is better)
"""
function fit_flux_two_step(flux, energies; kw...)
    # Use actual energy values as candidates (instrument channels)
    # Filter to reasonable range and ensure enough points on both sides
    sorted_energies = sort(unique(energies))
    Emins = filter(sorted_energies) do Emin
        n_low = sum(energies .<= Emin)
        n_high = sum(energies .> Emin)
        n_low >= 3 && n_high >= 5
    end

    isempty(Emins) && return nothing, nothing, nothing

    best_model = nothing
    best_score = Inf
    best_Emin = nothing
    best_flux_modeled = nothing

    for Emin in Emins
        try
            combined_model, flux_modeled = fit_flux_two_step(flux, energies, Emin; kw...)
            score = msd_log(flux, flux_modeled)
            if score < best_score
                best_score = score
                best_Emin = Emin
                best_flux_modeled = flux_modeled
                best_model = combined_model
            end
        catch
            continue
        end
    end
    return best_model, best_flux_modeled, best_Emin, best_score
end

export fit_two_flux

energies(x) = parent(x.dims[1].val)

function fit_two_flux(flux1, flux2; flux_threshold=100)
    flux = vcat(flux1, flux2)
    # Filter out NaN values and flux below threshold
    valid_idx = .!isnan.(flux) .& (flux .>= flux_threshold)
    clean_flux = flux[valid_idx]
    n_points = length(clean_flux)
    model, flux_modeled, Emin, score = fit_flux_two_step(parent(clean_flux), energies(clean_flux))
    success = !isnothing(model)
    return (; success, model, flux_modeled, n_points, Emin, score)
end
