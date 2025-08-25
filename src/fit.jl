export init_guess

include("model.jl")

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
    γ1, γ2 = -5.0, 3.0
    E0 = energies[i]
    Eb = 1.0e3
    A = flux[i] / exp(log_sbpl_model(E0, [1, γ1, γ2, Eb]))
    return [A, γ1, γ2, Eb]
end

function init_guess(::typeof(log_two_pop_model), energies, flux)
    i = argmax(flux)
    E0 = energies[i]
    Ec1 = E0
    Ec2 = 100.0e0
    γ1 = γ2 = 2.5
    A1 = flux[i] / E0^γ1 * exp(E0 / Ec1)
    A2 = A1 / 100
    return [A1, γ1, Ec1, A2, γ2, Ec2]
end

function init_guess(::typeof(log_plec_sbpl_model), energies, flux)
    p1 = [1.0e8, 3.0, 6.0e10]
    p2 = [8.0e7, -5.415023026784938, 2.9444485957997975, 63.69365352510846]
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

function fit_flux_two_step(flux, energies, Emin; kw...)
    alg = NonlinearSolve.TrustRegion()
    # first fit energy below 10^3 eV with log_plec_model
    m2 = log_sbpl_model_sciml

    flux_1 = flux[energies .< Emin]
    Es1 = energies[energies .< Emin]
    f1 = fit_powerlaw_exp(Es1, flux_1)

    # second fit the remaining flux of energy above 1 keV with log_sbpl_model
    δEs = energies[energies .>= Emin]
    δflux = flux[energies .>= Emin] - f1.(δEs)
    p2 = init_guess(m2, δEs, δflux)
    prob2 = NonlinearCurveFitProblem(m2, p2, δEs, log.(δflux))
    sol2 = solve(prob2; alg, kw...)

    f2 = SmoothBrokenPowerlaw(sol2.u..., 1.0)

    flux_1_modeled = f1.(Es1)
    flux_2_modeled = f2.(δEs) .+ f1.(δEs)

    return (Es1, f1), (δEs, f2), vcat(flux_1_modeled, flux_2_modeled)
end

"""
    _ssr(y, m)

Calculate the sum of squared residuals between y and m.
"""
_ssr(y, m) = sum(abs2, y .- m)

log_ssr(y, m) = sum(abs2, log.(y) .- log.(m))

# Method that optimizes Emin to minimize sum of squared residuals
function fit_flux_two_step(flux, energies; kw...)
    # Define the range of possible Emin values (exclude boundary energies)
    energy_range = energies[2:(end - 1)]  # Exclude first and last to ensure both regions have data

    best_ssr = Inf
    best_Emin = energy_range[1]
    best_result = nothing

    # Optimize Emin by trying different values and finding minimum SSR
    for Emin in energy_range
        # Ensure we have enough points in both regions
        # Need minimum points for fitting
        n_low = sum(energies .< Emin)
        n_high = sum(energies .>= Emin)
        if n_low < 3 || n_high < 4
            continue
        end

        try
            result = fit_flux_two_step(flux, energies, Emin; kw...)
            ssr = log_ssr(flux, result[3])

            # Update best if this is better
            if ssr < best_ssr
                best_ssr = ssr
                best_Emin = Emin
                best_result = result
            end
        catch
            continue
        end
    end

    return best_result, best_Emin, best_ssr
end

function fit_flux_two_step(flux::DimArray, args...; kw...)
    flux_clean = filter(!isnan, flux)
    fit_flux_two_step(parent(flux_clean), parent(flux_clean.dims[1].val), args...; kw...)
end