using NaNMath
import NaNMath as nm
using Statistics: mean
using LsqFit
export init_guess

# https://discourse.julialang.org/t/comparing-non-linear-least-squares-solvers/104752
# https://juliapackagecomparisons.github.io/comparisons/math/nonlinear_solvers/
# https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/index.html
# https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/plot_exp_cutoff_powerlaw.html
# https://fjebaker.github.io/SpectralFitting.jl/dev/

include("model.jl")

"""
    fit(PowerLawExpCutoff, E, y)

Fit the model f(E) = A * E^(-γ) * exp(-E/E_c) to data (E, y)
by minimizing ∑[ln(yᵢ) - ln f(Eᵢ)]².

Returns a NamedTuple with fields
- `A`: amplitude
- `γ`: power-law index
- `E_c`: cutoff energy
"""
function fit(::Type{<:PowerLawExpCutoff}, E, y)
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
    return PowerLawExpCutoff(A, γ, E_c)
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


# https://docs.sciml.ai/NonlinearSolve/stable/solvers/nonlinear_least_squares_solvers/
function sciml_log_fit(Model, E, y; kw...)
    f = (u, E) -> log_eval(Model(u...), E)
    p0 = init_guess(Model, E, y)
    prob = NonlinearCurveFitProblem(f, p0, E, nm.log.(y))
    sol = solve(prob; kw...)
    return Model(sol.u...)
end

# https://github.com/JuliaNLSolvers/LsqFit.jl
function LsqFit_log_fit(Model, E, y; kw...)
    f = (E, p) -> (m = Model(p...); log_eval(m, E))
    p0 = init_guess(Model, E, y)
    fit = LsqFit.curve_fit(f, E, nm.log.(y), p0; kw...)
    return Model(fit.param...)
end

function fit(::Type{<:SmoothBrokenPowerlaw}, E, y; kw...)
    f = log_sbpl_model_sciml
    alg = NonlinearSolve.TrustRegion()
    sol = sciml_log_fit(f, E, y; alg, kw...)
    return SmoothBrokenPowerlaw(sol.u..., 1.0)
end

function fit(M::Type{<:KappaDistribution}, E, y; kw...)
    # alg = NonlinearSolve.TrustRegion()
    # return sciml_log_fit(M, E, y; alg, kw...)
    LsqFit_log_fit(M, E, y; kw...)
end

function fit(M, flux; kw...)
    return fit(M, energies(flux), parent(flux); kw...)
end

function log_sbpl_model(E, p)
    A, γ1, γ2, Eb = p
    m = 1
    return log_sbpl.(E, A, γ1, γ2, Eb, m)
end

init_guess(M, flux) = init_guess(M, energies(flux), parent(flux))

function init_guess(::Type{<:KappaDistribution}, energies, flux)
    i = argmax(flux)
    κ0 = 5.0  # Typical kappa value
    E0 = energies[i]
    E_c0 = energies[i] / 1.1
    A = flux[i] / E0 * (1 + E0 / (κ0 * E_c0))^(κ0 + 1)
    return [A, κ0, E_c0]
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

for model in (:log_sbpl_model,)
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


function fit_flux_two_step(M1, M2, flux, energies, Emin; kw...)
    # first fit energy below Emin with model1 (default: PowerLawExpCutoff)
    flux_1 = flux[energies .<= Emin]
    Es1 = energies[energies .<= Emin]
    model1 = fit(M1, Es1, flux_1)

    # second fit the remaining flux of energy above Emin with model2 (default: SmoothBrokenPowerlaw)
    δEs = energies[energies .> Emin]
    δflux = flux[energies .> Emin] .- model1.(δEs)
    model2 = fit(M2, δEs, δflux)

    # Create combined model
    model = TwoStepModel(model1, model2, Emin)
    flux_modeled = model.(energies)

    return model, flux_modeled
end

"""
    msd_log(a, b)

Return the mean squared deviation between the log-transformed two arrays: `mean(abs2, log(a) - log(b))`.
"""
function msd_log(a, b)
    return mean(abs2, nm.log.(a) .- nm.log.(b))
end


"""
    fit(flux, energies; kw...)

Fit two-step model with optimized transition energy Emin.
Uses actual energy channel values as candidates for Emin.

Returns:
- TwoStepModel: Combined model structure
- best_Emin: Optimal transition energy
- best_score: Fit quality score (lower is better)
"""
function fit(t::Type{<:TwoStepModel{M1, M2}}, flux, energies; kw...) where {M1, M2}
    # Use actual energy values as candidates (instrument channels)
    # Filter to reasonable range and ensure enough points on both sides
    Emins = filter(energies) do Emin
        n_low = count(<=(Emin), energies)
        n_high = count(>(Emin), energies)
        Emin < 20 && n_low >= paramcount(M1) && n_high >= paramcount(M2)
    end
    isempty(Emins) && return nothing, nothing, nothing

    best_model = nothing
    best_score = Inf
    best_Emin = nothing
    best_flux_modeled = nothing

    for Emin in Emins
        try
            combined_model, flux_modeled = fit_flux_two_step(M1, M2, flux, energies, Emin; kw...)
            score = msd_log(flux, flux_modeled)
            if score < best_score
                best_score = score
                best_Emin = Emin
                best_flux_modeled = flux_modeled
                best_model = combined_model
            end
        catch
            # rethrow()
            continue
        end
    end
    return best_model, best_flux_modeled, best_Emin, best_score
end

export fit_two_flux

energies(x) = parent(x.dims[1].val)

function fit_two_flux(modelType, flux1, flux2; flux_threshold = 100)
    flux = vcat(flux1, flux2)
    # Filter out NaN values and flux below threshold
    valid_idx = @. !isnan(flux) & (flux >= flux_threshold)
    clean_flux = flux[valid_idx]
    n_points = length(clean_flux)
    model, flux_modeled, Emin, score = fit(modelType, parent(clean_flux), energies(clean_flux))
    success = !isnothing(model)
    return (; success, model, flux_modeled, n_points, Emin, score)
end

function fit_two_flux(flux1, flux2; flux_threshold = 100)
    modelType = TwoStepModel{PowerLawExpCutoff, SmoothBrokenPowerlaw}
    return fit_two_flux(modelType, flux1, flux2; flux_threshold)
end
