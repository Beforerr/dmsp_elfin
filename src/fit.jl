using NaNMath
import NaNMath as nm
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

# https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/plot_smooth_broken_powerlaw.html
struct SmoothBrokenPowerlaw{T}
    A::T
    γ1::T
    γ2::T
    Eb::T
    m::T
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

(f::SmoothBrokenPowerlaw)(E) = sbpl(E, f.A, f.γ1, f.γ2, f.Eb, f.m)

"""
    fit_powerlaw_exp(E, y)

Fit the model f(E) = A * E^(-γ) * exp(-E/E_c) to data (E, y)
by minimizing ∑[ln(yᵢ) - ln f(Eᵢ)]².

Returns a NamedTuple with fields
- `A`: amplitude
- `γ`: power-law index
- `E_c`: cutoff energy
"""
function fit_powerlaw_exp(E, y)
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

function log_plec_model(E, p)
    A, γ, Ec = p
    return nm.log(A) .- γ .* nm.log.(E) .- (E ./ Ec)
end

function log_plec_sbpl_model(E, p)
    A1, γ1, Ec1,    # PLEC params
    A2, γ2, γ3, Eb = p  # SBPL params
    m = 1

    # Component 1: PLEC
    f1 = @. A1 * E^(-γ1) * exp(-E / Ec1)

    # Component 2: SBPL
    x = E ./ Eb
    f2 = @. A2 * x^(-γ2) * (1 + x^m)^((γ2 - γ3) / m)

    return log10.(f1 .+ f2)
end

function log_two_pop_model(p, E)
    A1, γ1, Ec1, A2, γ2, Ec2 = p
    f1 = A1 .* E .^ (-γ1) .* exp.(-E ./ Ec1)
    f2 = A2 .* E .^ (-γ2) .* exp.(-E ./ Ec2)
    flux = f1 .+ f2
    return NaNMath.log10.(flux)
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
    γ1, γ2 = -5., 3.
    E0 = energies[i]
    Eb = 1e3
    A = flux[i] / exp(log_sbpl_model(E0, [1, γ1, γ2, Eb]))
    return [A, γ1, γ2, Eb]
end

function init_guess(::typeof(log_two_pop_model), energies, flux)
    i = argmax(flux)
    E0 = energies[i]
    Ec1 = E0
    Ec2 = 100E0
    γ1 = γ2 = 2.5
    A1 = flux[i] / E0^γ1 * exp(E0 / Ec1)
    A2 = A1 / 100
    return [A1, γ1, Ec1, A2, γ2, Ec2]
end

function init_guess(::typeof(log_plec_sbpl_model), energies, flux)
    p1 = [1e8, 3., 6e10]
    p2 = [8e7, -5.415023026784938, 2.9444485957997975, 63.69365352510846]
    return [p1..., p2...]
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


function fit_flux_two_step(flux, energies, Emin=1; kw...)
    alg = NonlinearSolve.TrustRegion()
    # first fit energy below 10^3 eV with log_plec_model
    m2 = log_sbpl_model_sciml

    flux_1 = flux[energies.<Emin]
    Es1 = energies[energies.<Emin]
    f1 = fit_powerlaw_exp(Es1, flux_1)

    # second fit the remaining flux of energy above 1 keV with log_sbpl_model
    δEs = energies[energies.>=Emin]
    δflux = flux[energies.>=Emin] - f1.(δEs)
    p2 = init_guess(m2, δEs, δflux)
    prob2 = NonlinearCurveFitProblem(m2, p2, δEs, log.(δflux))
    sol2 = solve(prob2; alg, kw...)

    f2 = SmoothBrokenPowerlaw(sol2.u..., 1.)

    flux_1_modeled = f1.(Es1)
    flux_2_modeled = f2.(δEs) .+ f1.(δEs)

    (Es1, f1), (δEs, f2), vcat(flux_1_modeled, flux_2_modeled)
end
