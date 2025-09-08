using NaNMath
import NaNMath as nm
using Statistics: mean
using LsqFit
using StaticArrays
export init_guess

# https://discourse.julialang.org/t/comparing-non-linear-least-squares-solvers/104752
# https://juliapackagecomparisons.github.io/comparisons/math/nonlinear_solvers/
# https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/index.html
# https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/plot_exp_cutoff_powerlaw.html
# https://fjebaker.github.io/SpectralFitting.jl/dev/

include("model.jl")
include("JSO.jl")
include("SciML.jl")

# https://github.com/JuliaNLSolvers/LsqFit.jl
function LsqFit_log_fit(Model, E, y; kw...)
    f = (E, p) -> (m = Model(p); log_eval.(m, E))
    p0 = init_guess(Model, E, y)
    log_y = log.(y)
    try
        fit = LsqFit.curve_fit(f, E, log_y, p0; kw...)
        return Model(fit.param)
    catch e
        @info log_y p0
        throw(e)
    end
end

init_guess(M, flux) = init_guess(M, energies(flux), parent(flux))

function unbound_fit(M::Type{<:PowerLawExpCutoff}, E, y)
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
    return M(A, γ, E_c)
end

"""
    init_guess(PowerLawExpCutoff, E, y)

Note that the result may not be a valid PowerLawExpCutoff model (i.e. E_c <= 0)
"""
function init_guess(M::Type{<:PowerLawExpCutoff}, E, y)
    m = unbound_fit(M, E, y)
    return if m.E_c > 0
        raw_vec(m)
    else
        _init_guess(M, E, y)
    end
end

function init_guess(::Type{<:PowerLawExpCutoff2}, E, y)
    m = unbound_fit(PowerLawExpCutoff, E, y)
    A, γ, E_c = if m.E_c > 0
        raw_vec(m)
    else
        _init_guess(PowerLawExpCutoff, E, y)
    end
    return [log(A), γ, log(E_c)]
end

function _init_guess(::Type{<:PowerLawExpCutoff}, E, y)
    i = argmax(y)
    γ = -2
    E_c = E[i]
    A = y[i] * E[i]^γ * exp(E[i] / E_c)
    return [A, γ, E_c]
end

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

"""
    fit(PowerLawExpCutoff, E, y) => PowerLawExpCutoff(A, γ, E_c)

Fit the model f(E) = A * E^(-γ) * exp(-E/E_c) to data (E, y)
by minimizing ∑[ln(yᵢ) - ln f(Eᵢ)]².
"""
function fit(M::Type{<:PowerLawExpCutoff}, E, y)
    # return M(init_guess(M, E, y))
    lower = SA[0, -Inf, 0]
    method = :sciml
    if method == :sciml
        alg = NonlinearSolve.TrustRegion()
        return sciml_log_fit(M, E, y; alg)
    elseif method == :lsqfit
        return LsqFit_log_fit(M, E, y; lower)
    end
end

"""
    fit(PowerLaw, E, y) => f::PowerLaw(A, γ)

Fit the model f(E) = A * E^(-γ) to data (E, y) by minimizing ∑[ln(yᵢ) - ln f(Eᵢ)]².
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

function fit(M::Type{<:KappaDistribution}, E, y; method = :sciml, kw...)
    return if method == :sciml
        alg = NonlinearSolve.TrustRegion() # this is much faster than the default `FastShortcutNLLSPolyalg()` and `LeastSquaresOptimJL()`
        sciml_log_fit(M, E, y; alg, kw...)
    elseif method == :jso
        jso_nls_fit(M, E, y; kw...)
    elseif method == :lsqfit
        LsqFit_log_fit(M, E, y; kw...)
    end
end

function fit(M, flux; kw...)
    return fit(M, energies(flux), parent(flux); kw...)
end

export remove_nan, fit_flux_two_step, jso_nls_fit

# TODO: ignore the last channel of dmsp_flux

"""
    fit_flux_two_step(M1, M2, energies, flux, E_tran; kw...)

First fit high energy part (E > E_tran) with M2
Then fit low energy part (E <= E_tran) with M1
"""
function fit_flux_two_step(M1, M2, energies, flux, E_tran; kw...)
    high_energy_idx = energies .> E_tran
    flux_high = flux[high_energy_idx]
    Es_high = energies[high_energy_idx]
    model2 = fit(M2, Es_high, flux_high; kw...)

    # fit low energy part
    low_energy_idx = energies .<= E_tran
    flux_low = flux[low_energy_idx]
    Es_low = energies[low_energy_idx]

    alg = NonlinearSolve.TrustRegion()
    u1_0 = init_guess(M1, Es_low, flux_low)
    function f(u, E)
        m1 = M1(u)
        return @. nm.log(m1(E) + model2(E))
    end
    prob = NonlinearCurveFitProblem(f, u1_0, energies, log.(flux))
    sol = solve(prob, alg; kw...)

    model1 = M1(sol.u)
    model = TwoStepModel(model1, model2, E_tran)
    score = mean(abs2, sol.resid)
    return model, score
end


# This does not work very well as quite nonlinear without constraints
function fit_flux_one_step(M1, M2, energies, flux, E_tran; kw...)
    high_energy_idx = energies .> E_tran
    flux_high = flux[high_energy_idx]
    Es_high = energies[high_energy_idx]
    model2 = fit(M2, Es_high, flux_high; kw...)

    # fit low energy part
    low_energy_idx = energies .<= E_tran
    flux_low = flux[low_energy_idx]
    Es_low = energies[low_energy_idx]
    δflux_low = flux_low .- model2.(Es_low)
    if any(<(0), δflux_low) # Not a valid fit for low energy part
        return missing, Inf
    end

    alg = NonlinearSolve.TrustRegion()
    u1_0 = init_guess(M1, Es_low, δflux_low)
    f = (u, E) -> begin
        m1 = M1(u[1:paramcount(M1)])
        m2 = M2(u[(paramcount(M1) + 1):end])
        @. nm.log(m1(E) + m2(E))
    end
    u0 = vcat(u1_0, raw_vec(model2))
    prob = NonlinearCurveFitProblem(f, u0, energies, log.(flux))
    sol = solve(prob, alg; kw...)
    model = TwoStepModel(M1(sol.u[1:paramcount(M1)]), M2(sol.u[(paramcount(M1) + 1):end]), E_tran)
    score = msd_log(model, energies, flux)
    return model, score
end

"""
    msd_log(a, b)

Return the mean squared deviation between the log-transformed two arrays: `mean(abs2, log(a) - log(b))`.
"""
msd_log(a, b) = mean(abs2, log.(a) .- log.(b))
msd_log(model, x, y) = mean(abs2, log.(y) .- log.(model.(x)))


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
    best_model = missing
    best_score = Inf

    n_low = paramcount(M1)
    n_high = paramcount(M2)

    for E_tran in energies
        if E_tran > 20 || count(<=(E_tran), energies) < n_low || count(>(E_tran), energies) < n_high
            continue
        end
        combined_model, score = fit_flux_two_step(M1, M2, energies, flux, E_tran; kw...)
        verbose && @info combined_model
        if score < best_score
            best_score = score
            best_model = combined_model
        end
    end
    return best_model, best_score
end

export fit_two_flux

function fit_two_flux(modelType, flux1, flux2; flux_threshold = 200, kw...)
    flux = vcat(flux1, flux2)

    # Filter out NaN values and flux below threshold
    valid_idx = @. !isnan(flux) & (flux >= flux_threshold)
    clean_flux = flux[valid_idx]
    n_points = length(clean_flux)
    model, score = fit(modelType, parent(clean_flux), energies(clean_flux); kw...)
    return (; success = !ismissing(model), model, n_points, score)
end

function fit_two_flux(flux1, flux2; kw...)
    modelType = TwoStepModel{PowerLawExpCutoff2, KappaDistribution}
    return fit_two_flux(modelType, flux1, flux2; kw...)
end
