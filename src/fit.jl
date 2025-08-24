using NaNMath
import NaNMath as nm
export init_guess

include("model.jl")

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

export fit_row_parameters

energies(x) = parent(x.dims[1].val)

function fit_row_parameters(flux1, flux2; mlat = nothing)
    ff = vcat(parent(flux1), parent(flux2))
    ee = vcat(energies(flux1), energies(flux2))
    ff, ee = remove_nan(ff, ee)
    n_points = length(ff)

    Emin = 90

    if n_points < 5  # Need minimum points for fitting
        return (; success = false, flux_modeled = nothing, params = nothing, n_points)
    end
    try
        (Es1, f1), (δEs, f2), flux_modeled = fit_flux_two_step(ff, ee, Emin)
        return (; success = true, flux_modeled, params = (f1, f2, Emin), n_points)
    catch e
        @warn "Failed to fit MLAT $(mlat): $e"
        return (; success = false, flux_modeled = nothing, params = nothing, n_points)
    end
end
