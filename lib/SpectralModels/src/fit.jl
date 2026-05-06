# https://discourse.julialang.org/t/comparing-non-linear-least-squares-solvers/104752
# https://juliapackagecomparisons.github.io/comparisons/math/nonlinear_solvers/
# https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/index.html

using Statistics: mean
export init_guess, fit_two_step, fit_two_flux
using SciMLBase: FullSpecialize
using NonlinearSolveFirstOrder: TrustRegion

# IIP + FullSpecialize avoids FunctionWrappersWrapper (which locks IIP to chunk_size=1)
const _FS = FullSpecialize

energies(x) = hasfield(typeof(x), :dims) ? Base.parent(x.dims[1].val) : x

init_guess(M, flux) = init_guess(M, energies(flux), Base.parent(flux))

function init_guess(M::Type{<:AbstractPowerLawExpCutoff}, E, y)
    function _init_guess(E, y)
        i = argmax(y)
        γ = -1.5 # Like Maxwellian distribution
        E_c = E[i]
        A = y[i] * E[i]^γ * exp(E[i] / E_c)
        return PowerLawExpCutoff(A, γ, E_c)
    end

    m = unbound_fit(PowerLawExpCutoff, E, y)
    return M(m.E_c > 0 ? m : _init_guess(E, y))
end

function init_guess(M::Type{<:AbstractKappaDistribution}, energies, flux)
    function _init_guess(energies, flux)
        i = argmax(flux)
        κ0 = 5.0  # Typical kappa value
        E0 = energies[i]
        E_c0 = energies[i] / 1.1
        A = flux[i] / E0 * (1.0 + E0 / (κ0 * E_c0))^(κ0 + 1.0)
        return KappaDistribution(A, κ0, E_c0)
    end

    return M(_init_guess(energies, flux))
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
    return if method == :sciml
        alg = TrustRegion()
        sciml_log_fit(M, E, y; alg)
    elseif method == :lsqfit
        LsqFit_log_fit(M, E, y; lower)
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

function fit(M::Type{<:AbstractKappaDistribution}, E, y; method = :sciml, kw...)
    return if method == :sciml
        alg = TrustRegion() # this is much faster than the default `FastShortcutNLLSPolyalg()` and `LeastSquaresOptimJL()`
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

# convert model to avoid repeated transformation
_eval_model(m) = m
_eval_model(m::PowerLawExpCutoff2) = PowerLawExpCutoff(A(m), m.γ, E_c(m))
_eval_model(m::TransformKappaDistribution) = KappaDistribution(m)


"""
    fit_two_step(m1_init, m2_init, x, y; high, kw...)

Fit M2 (seeded from `m2_init`) to `x[high], y[high]`, then fit M1 (seeded from `m1_init`)
to all of (x, y) with M2 as additive background. Returns `(model1, model2, score)`.
"""
function fit_two_step(m1_init::M1, m2_init::M2, x, y; high, logy = log.(y), refine = false, kw...) where {M1, M2}
    alg = TrustRegion()

    model2 = fit(M2, x[high], y[high]; u0 = raw_vec(m2_init), kw...)
    m2_vals = model2.(x)

    model1 = begin
        f! = (r, u, E) -> (m1 = _eval_model(M1(u)); @. r = log(m1(E) + m2_vals) - logy)
        nf = NonlinearFunction{true, _FS}(f!; resid_prototype = similar(x))
        prob = NonlinearLeastSquaresProblem(nf, raw_vec(m1_init), x)
        sol = solve(prob, alg; kw...)
        M1(sol.u)
    end

    if refine
        n1 = paramcount(M1)
        u0_joint = vcat(raw_vec(model1), raw_vec(model2))
        f_joint! = (r, u, E) -> (
            m1 = _eval_model(M1(@view u[1:n1]));
            m2 = _eval_model(M2(@view u[(n1 + 1):end]));
            @. r = log(m1(E) + m2(E)) - logy
        )
        nf_joint = NonlinearFunction{true, _FS}(f_joint!; resid_prototype = similar(x))
        prob_joint = NonlinearLeastSquaresProblem(nf_joint, u0_joint, x)
        sol = solve(prob_joint, alg; kw...)
        model1 = M1(@view sol.u[1:n1])
        model2 = M2(@view sol.u[(n1 + 1):end])
    end

    score = mean(abs2, sol.resid)
    return model1, model2, score
end

"""
    fit(TwoStepModel, flux, energies; kw...)

Fit two-step model with optimized transition energy.

Returns:
- TwoStepModel: Combined model structure
- best_score: Fit quality score (lower is better)
"""
function fit(t::Type{<:TwoStepModel{M1, M2}}, flux, energies; EMAX = 60, kw...) where {M1, M2}
    if !issorted(energies)
        p = sortperm(energies)
        return fit(t, flux[p], energies[p]; EMAX, kw...)
    end

    best_model, best_score = missing, Inf

    n_low = paramcount(M1)
    n_high = paramcount(M2)
    n = length(energies)
    logy = log.(flux)

    # Warm-start: store fitted models and reuse as init for adjacent E_tran candidates
    m1_warm = nothing
    m2_warm = nothing
    @views for (i, E_tran) in enumerate(energies)
        if E_tran > EMAX || i < n_low || n - i < n_high
            continue
        end
        high = (i + 1):n
        m1_init = @something m1_warm init_guess(M1, energies[1:i], flux[1:i])
        m2_init = @something m2_warm init_guess(M2, energies[high], flux[high])
        model1, model2, score = fit_two_step(m1_init, m2_init, energies, flux; high, logy, kw...)
        if score < best_score
            best_model, best_score = TwoStepModel(model1, model2, E_tran), score
        end
        m1_warm = model1
        m2_warm = model2
    end
    return best_model, best_score
end

const FLUX_THRESHOLD = 150

function fit_two_flux(modelType, flux1, flux2; threshold = FLUX_THRESHOLD, kw...)
    flux = vcat(flux1, flux2)

    # Filter out NaN values and flux below threshold
    valid_idx = @. !isnan(flux) & (flux >= threshold)
    clean_flux = flux[valid_idx]
    n_points = length(clean_flux)
    model, score = fit(modelType, Base.parent(clean_flux), energies(clean_flux); kw...)
    return (; success = !ismissing(model), model, n_points, score)
end

function fit_two_flux(flux1, flux2; modelType = TwoStepModel{PowerLawExpCutoff2, TransformKappaDistribution}, kw...)
    return fit_two_flux(modelType, flux1, flux2; kw...)
end
