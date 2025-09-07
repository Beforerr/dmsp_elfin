# https://fjebaker.github.io/SpectralFitting.jl/dev/models/using-models/#SpectralFitting.AbstractSpectralModel
import Base: iterate, show

export SpectralModel, TwoStepModel, PowerLaw, PowerLawExpCutoff, SmoothBrokenPowerlaw, KappaDistribution

abstract type SpectralModel{T} end

paramcount(m::Type{<:SpectralModel}) = fieldcount(m)
paramcount(m::SpectralModel) = paramcount(m)

function Base.iterate(m::T, state = 1) where {T <: SpectralModel}
    return state > fieldcount(T) ? nothing : (getfield(m, state), state + 1)
end

Base.broadcastable(o::SpectralModel) = Ref(o)

function Base.show(io::IO, m::T) where {T <: SpectralModel}
    print(io, Base.typename(T).name, "(")
    for fn in fieldnames(T)
        fv = getfield(m, fn)
        print(io, fn, "=", @sprintf("%.2g", fv))
        fn != last(fieldnames(T)) && print(io, ", ")
    end
    print(io, ")")
    return
end

"""
    T(u::AbstractVector) where {T <: SpectralModel}

Construct a spectral model from a parameter vector without splatting for performance.
"""
@generated function (::Type{T})(u::AbstractVector) where {T <: SpectralModel}
    n = fieldcount(T)
    args = [:(u[$i]) for i in 1:n]
    return quote
        @assert length(u) == $n
        @inbounds $(T)($(args...))
    end
end

"""
    raw_vec(m)

Extract raw parameter vector from model object.
"""
raw_vec(m::SpectralModel) = [getfield(m, fn) for fn in fieldnames(typeof(m))]

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
@kwdef struct PowerLawExpCutoff{T} <: SpectralModel{T}
    A::T
    γ::T
    E_c::T
end

(m::PowerLawExpCutoff)(E) = m.A * E^(-m.γ) * exp(-E / m.E_c)

log_eval(m::PowerLawExpCutoff, E) = nm.log(m.A) - m.γ * log(E) - E / m.E_c

"""
Kappa distribution spectral model.

```math
f(E) = A * E * (1 + E/(κ*E_c))^(-κ-1)
```

Where κ is the kappa parameter controlling the suprathermal tail.
"""
@kwdef struct KappaDistribution{T} <: SpectralModel{T}
    A::T
    κ::T
    E_c::T
end

(m::KappaDistribution)(E) = m.A * E * (1 + E / (m.κ * m.E_c))^(-m.κ - 1)

log_eval(m::KappaDistribution, E) = nm.log(m.A) + log(E) + nm.log(1 + E / (m.κ * m.E_c)) * (-m.κ - 1)

function log_jacobian(m::KappaDistribution, E)
    A = m.A
    κ = m.κ
    E_c = m.E_c
    u = E / (κ * E_c)           # u = E/(κ E_c)
    g = nm.log(1 + u)

    dA = 1 / A
    dκ = -g + (κ + 1) * E / (κ * (κ * E_c + E))
    dEc = (κ + 1) * E / (E_c * (κ * E_c + E))

    return [dA, dκ, dEc]
end

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
    return nm.log(A) - γ1 * nm.log(E) + ((γ1 - γ2) / m) * nm.log(1 + x^m)
end

(m::SmoothBrokenPowerlaw)(E) = exp(log_sbpl(E, m...))

# Fixed m = 1
struct SmoothBrokenPowerlawFixed{T} <: SpectralModel{T}
    A::T
    γ1::T
    γ2::T
    Eb::T
end


(m::SmoothBrokenPowerlawFixed)(E) = exp(log_sbpl(E, m..., 1))

"""
    TwoStepModel{M1,M2,T}

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
# With PowerLawExpCutoff and SmoothBrokenPowerlaw
plec = PowerLawExpCutoff(A, γ, E_c)
sbpl = SmoothBrokenPowerlaw(A, γ1, γ2, Eb, m)
model = TwoStepModel(plec, sbpl, 100.0)

# With any callable models
gaussian = x -> exp(-x^2)
exponential = x -> exp(-x)
model = TwoStepModel(gaussian, exponential, 2.0)
```
"""
struct TwoStepModel{M1, M2, T}
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

function (m::TwoStepModel)(E)
    return m.model2(E) + m.model1(E)
end
