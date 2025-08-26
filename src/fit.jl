using NaNMath
import NaNMath as nm
using Statistics: mean
using LsqFit
using JSOSolvers, CaNNOLeS
using ADNLPModels
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
    f = (u, E) -> (m = Model(u...); log_eval.(m, E))
    p0 = init_guess(Model, E, y)
    prob = NonlinearCurveFitProblem(f, p0, E, nm.log.(y))
    sol = solve(prob; kw...)
    return Model(sol.u...)
end

# https://github.com/JuliaNLSolvers/LsqFit.jl
function LsqFit_log_fit(Model, E, y; kw...)
    f = (E, p) -> (m = Model(p...); log_eval.(m, E))
    p0 = init_guess(Model, E, y)
    fit = LsqFit.curve_fit(f, E, nm.log.(y), p0; kw...)
    return Model(fit.param...)
end

# JSO (JuliaSmoothOptimizers) interface for nonlinear least squares
"""
    jso_nls_fit(Model, E, y; solver=cannoles, bounds=nothing, kw...)

Fit spectral model using JSO nonlinear least squares solvers.

Uses ADNLPModels.jl for problem modeling and CaNNOLeS.jl for solving.

# Example
```julia
using CaNNOLeS
model = jso_nls_fit(KappaDistribution, energies, flux)

# Or with JSOSolvers
using JSOSolvers
model = jso_nls_fit(KappaDistribution, energies, flux, solver=trunk)
```

# References
- https://jso.dev/ADNLPModels.jl/dev/
- https://jso.dev/JSOSolvers.jl/latest/solvers/
"""
function jso_nls_fit(Model, E, y; bounds = nothing, kw...)
    # Residual function for nonlinear least squares
    observed = nm.log.(y)
    function F(x)
        m = Model(x...)
        return observed .- log_eval.(m, E)
    end

    p0 = init_guess(Model, E, y)

    # Create constrained ADNLPModel
    nls = if bounds !== nothing
        lvar, uvar = bounds
        ADNLSModel(F, p0, length(E), lvar, uvar)
    else
        ADNLSModel(F, p0, length(E))
    end

    # Solve using JSO solver
    stats = if isnothing(bounds)
        trunk(nls; kw...)
    else
        tron(nls)
    end

    return Model(stats.solution...)
end

function fit(::Type{<:SmoothBrokenPowerlawFixed}, E, y; kw...)
    f = log_sbpl_model_sciml
    alg = NonlinearSolve.TrustRegion()
    sol = sciml_log_fit(f, E, y; alg, kw...)
    return SmoothBrokenPowerlawFixed(sol.u...)
end

function fit(M::Type{<:KappaDistribution}, E, y; method = :sciml, kw...)
    return if method == :jso
        jso_nls_fit(M, E, y; kw...)
    elseif method == :sciml
        alg = NonlinearSolve.TrustRegion()
        sciml_log_fit(M, E, y; alg, kw...)
    else  # Default: LsqFit
        LsqFit_log_fit(M, E, y; kw...)
    end
end

function fit(M, flux; kw...)
    return fit(M, energies(flux), parent(flux); kw...)
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

function init_guess(::Type{<:SmoothBrokenPowerlawFixed}, energies, flux)
    i = argmax(flux)
    γ1, γ2 = -5.0, 3.0
    E0 = energies[i]
    Eb = 1.0e3
    A = flux[i] / exp(log_sbpl_model(E0, [1, γ1, γ2, Eb]))
    return [A, γ1, γ2, Eb]
end

using LeastSquaresOptim, NonlinearSolve
using CurveFit
using ADNLPModels, CaNNOLeS, JSOSolvers

export remove_nan, fit_flux_two_step, jso_nls_fit, jso_constrained_nls_fit, jso_trunk_fit

# TODO: ignore the last channel of dmsp_flux
"""
    fit_flux_two_step(M1, M2, flux, energies, Emin; kw...)

First fit high energy part (E > Emin) with model2
Then fit low energy (cold) part after subtracting high energy contribution
"""
function fit_flux_two_step(M1, M2, flux, energies, Emin; kw...)
    high_energy_idx = energies .> Emin
    flux_high = flux[high_energy_idx]
    Es_high = energies[high_energy_idx]
    model2 = fit(M2, Es_high, flux_high; kw...)

    low_energy_idx = energies .<= Emin
    flux_low = flux[low_energy_idx]
    Es_low = energies[low_energy_idx]
    δflux_low = flux_low .- model2.(Es_low)
    model1 = fit(M1, Es_low, δflux_low)

    model = TwoStepModel(model1, model2, Emin)
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
    fit(flux, energies; kw...)

Fit two-step model with optimized transition energy Emin.
Uses actual energy channel values as candidates for Emin.

Returns:
- TwoStepModel: Combined model structure
- best_Emin: Optimal transition energy
- best_score: Fit quality score (lower is better)
"""
function fit(t::Type{<:TwoStepModel{M1, M2}}, flux, energies; verbose = false, kw...) where {M1, M2}
    # Use actual energy values as candidates (instrument channels)
    # Filter to reasonable range and ensure enough points on both sides
    best_model = nothing
    best_score = Inf
    best_Emin = nothing
    best_flux_modeled = nothing

    n_low = paramcount(M1)
    n_high = paramcount(M2)

    for Emin in energies
        if Emin > 20 || count(<=(Emin), energies) < n_low || count(>(Emin), energies) < n_high
            continue
        end
        try
            combined_model, flux_modeled = fit_flux_two_step(M1, M2, flux, energies, Emin; kw...)
            score = msd_log(flux, flux_modeled)
            verbose && @info combined_model
            if score < best_score
                best_score = score
                best_Emin = Emin
                best_flux_modeled = flux_modeled
                best_model = combined_model
            end
        catch e
            verbose && @info "Failed to fit Emin = $Emin" e
            continue
        end
    end
    return best_model, best_flux_modeled, best_Emin, best_score
end

export fit_two_flux

function fit_two_flux(modelType, flux1, flux2; flux_threshold = 100, kw...)
    flux = vcat(flux1, flux2)
    # Filter out NaN values and flux below threshold
    valid_idx = @. !isnan(flux) & (flux >= flux_threshold)
    clean_flux = flux[valid_idx]
    n_points = length(clean_flux)
    model, flux_modeled, Emin, score = fit(modelType, parent(clean_flux), energies(clean_flux); kw...)
    success = !isnothing(model)
    return (; success, model, flux_modeled, n_points, Emin, score)
end

function fit_two_flux(flux1, flux2; kw...)
    modelType = TwoStepModel{PowerLawExpCutoff, KappaDistribution}
    return fit_two_flux(modelType, flux1, flux2; kw...)
end
