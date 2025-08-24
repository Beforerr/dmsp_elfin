using NaNMath
import NaNMath as nm
using Statistics: mean
export init_guess, TwoStepModel

# https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/plot_exp_cutoff_powerlaw.html

"""
Power-law model with exponential cutoff.

```math
f(E) = A * E^(-γ) * exp(-E/E_c)
```
"""
struct PowerLawExp{T}
    A::T
    γ::T
    E_c::T
end

(f::PowerLawExp)(E) = f.A * E^(-f.γ) * exp(-E / f.E_c)

"""
    TwoStepModel{T,M1,M2}

General two-step combined model with a transition energy Emin.

The combined model is defined as:
- For E ≤ Emin: f(E) = model1(E)
- For E > Emin: f(E) = model1(E) + model2(E)

# Fields
- `model1::M1`: First model (typically for low energy region)
- `model2::M2`: Second model (typically for high energy region)
- `Emin::T`: Transition energy between the two components

# Usage
```julia
model = TwoStepModel(model1, model2, Emin)
flux = model(energy)  # Evaluate at any energy
```



# Examples
```julia
# With PowerLawExp and SmoothBrokenPowerlaw
plec = PowerLawExp(A, γ, E_c)
sbpl = SmoothBrokenPowerlaw(A, γ1, γ2, Eb, m)
model = TwoStepModel(plec, sbpl, 100.0)

# With any callable models
gaussian = x -> exp(-x^2)
exponential = x -> exp(-x)
model = TwoStepModel(gaussian, exponential, 2.0)
```
"""
struct TwoStepModel{T, M1, M2}
    model1::M1
    model2::M2
    Emin::T
end

# Indexing interface for backward compatibility
Base.getindex(model::TwoStepModel, i::Integer) = if i == 1
    model.model1
elseif i == 2
    model.model2
elseif i == 3
    model.Emin
else
    throw(ArgumentError("Index out of bounds"))
end

function (model::TwoStepModel)(E)
    model1_flux = model.model1(E)
    if E <= model.Emin
        return model1_flux
    else
        model2_flux = model.model2(E)
        return model1_flux + model2_flux
    end
end

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
"""
    SmoothBrokenPowerlaw{T}

Smooth broken power-law model. 

```math
f(E) = A * E^(-γ1) * (1 + (E/Eb)^m)^((γ2 - γ1)/m)
```
"""
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
    flux_1 = flux[energies .<= Emin]
    Es1 = energies[energies .<= Emin]
    f1 = fit(PowerLawExp, Es1, flux_1)

    # second fit the remaining flux of energy above Emin with `SmoothBrokenPowerlaw`
    δEs = energies[energies .> Emin]
    δflux = flux[energies .> Emin] - f1.(δEs)
    f2 = fit(SmoothBrokenPowerlaw, δEs, δflux)

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
        n_low >= 3 && n_high >= 3
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

export fit_row_parameters

energies(x) = parent(x.dims[1].val)

function fit_row_parameters(flux1, flux2)
    ff = vcat(parent(flux1), parent(flux2))
    ee = vcat(energies(flux1), energies(flux2))
    ff, ee = remove_nan(ff, ee)
    n_points = length(ff)
    model, flux_modeled, Emin, score = fit_flux_two_step(ff, ee)
    success = !isnothing(model)
    return (; success, model, flux_modeled, n_points, Emin, score)
end
