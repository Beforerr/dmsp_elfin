# https://fjebaker.github.io/SpectralFitting.jl/dev/models/using-models/#SpectralFitting.AbstractSpectralModel
import Base: iterate

abstract type SpectralModel{T} end

function Base.iterate(m::T, state = 1) where {T <: SpectralModel}
    return state > fieldcount(T) ? nothing : (getfield(m, state), state + 1)
end

"""
Power-law model.

```math
f(E) = A * E^(-γ)
```
"""
struct PowerLaw{T} <: SpectralModel{T}
    A::T
    γ::T
end

(f::PowerLaw)(E) = f.A * E^(-f.γ)

"""
Power-law model with exponential cutoff.

```math
f(E) = A * E^(-γ) * exp(-E/E_c)
```
"""
struct PowerLawExp{T} <: SpectralModel{T}
    A::T
    γ::T
    E_c::T
end

(f::PowerLawExp)(E) = f.A * E^(-f.γ) * exp(-E / f.E_c)


"""
    SmoothBrokenPowerlaw{T}

Smooth broken power-law model.

```math
f(E) = A * E^(-γ1) * (1 + (E/Eb)^m)^((γ2 - γ1)/m)
```

# References
- [gammapy](https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/plot_smooth_broken_powerlaw.html)
"""
struct SmoothBrokenPowerlaw{T} <: SpectralModel{T}
    A::T
    γ1::T
    γ2::T
    Eb::T
    m::T
end

function log_sbpl(E, A, γ1, γ2, Eb, m)
    x = E / Eb
    return log(A) - γ1 * log(E) + ((γ1 - γ2) / m) * log(1 + x^m)
end

(m::SmoothBrokenPowerlaw)(E) = exp(log_sbpl(E, m...))


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
