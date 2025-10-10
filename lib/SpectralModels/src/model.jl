"""
    SpectralModel{T}

Abstract base type for all spectral models.

Features:
- Callable: `model(E)` evaluates flux at energy E
- Parameter extraction: `raw_vec(model)`
"""
abstract type SpectralModel{T} end

"""
    paramcount(model)
    paramcount(ModelType)

Return the number of parameters in the spectral model.
"""
paramcount(m::Type{<:SpectralModel}) = fieldcount(m)
paramcount(m::SpectralModel) = paramcount(typeof(m))

"""
    n_flux(model, Emin, Emax)
Compute number flux (particle count) by integrating the spectral model from Emin to Emax.

# Examples
```julia
model = PowerLawExpCutoff2(logA=10.0, γ=2.0, logE_c=3.0)
flux = n_flux(model, 10.0, 100.0)  # Total particles between 10-100 keV
```
"""
n_flux(m, Emin, Emax) = _integral0(m, Emax) - _integral0(m, Emin)
n_flux(Emin, Emax) = m -> n_flux(m, Emin, Emax)

"""
    e_flux(model, Emin, Emax) -> Float64

Compute energy flux by integrating E·f(E) from Emin to Emax.

# Examples
```julia
model = PowerLawExpCutoff2(logA=10.0, γ=2.0, logE_c=3.0)
energy_flux = e_flux(model, 10.0, 100.0)  # Total energy between 10-100 keV
```
"""
e_flux(m, Emin, Emax) = _integral1(m, Emax) - _integral1(m, Emin)
e_flux(Emin, Emax) = m -> e_flux(m, Emin, Emax)

# Make models iterable over their parameters
function Base.iterate(m::T, state = 1) where {T <: SpectralModel}
    return state > fieldcount(T) ? nothing : (getfield(m, state), state + 1)
end

# Make models broadcastable
Base.broadcastable(o::SpectralModel) = Ref(o)

# Pretty printing
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
    (::Type{T})(u::AbstractVector) where {T <: SpectralModel}

Construct a spectral model from a parameter vector.

This method avoids splatting for better performance during optimization.
"""
@generated function (::Type{T})(u::AbstractVector) where {T <: SpectralModel}
    n = fieldcount(T)
    args = [:(u[$i]) for i in 1:n]
    return quote
        @assert length(u) == $n "Expected $($n) parameters, got $(length(u))"
        @inbounds $(T)($(args...))
    end
end

"""
    raw_vec(model::SpectralModel)

Extract parameter vector from model object.

# Examples
```julia
model = PowerLaw(1e5, 2.0)
params = raw_vec(model)  # [1e5, 2.0]
```
"""
raw_vec(m::SpectralModel) = [getfield(m, fn) for fn in fieldnames(typeof(m))]

"""
    math_show(model)

Return LaTeX representation of the model with parameter values.
"""
math_show(m) = string(m)


# ============================================================================
# Broken Power-Law Models
# ============================================================================

"""
    SmoothBrokenPowerlaw{T} <: SpectralModel{T}

Smooth broken power-law model with variable smoothness.

# Model
```math
f(E) = A ⋅ E^{-γ_1} ⋅ [1 + (E/E_b)^m]^{(γ_2 - γ_1)/m}
```

# Fields
- `A::T`: Normalization amplitude
- `γ1::T`: Low-energy spectral index
- `γ2::T`: High-energy spectral index
- `Eb::T`: Break energy
- `m::T`: Smoothness parameter (m=1 gives smooth transition)

# References
[Gammapy smooth broken powerlaw](https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/plot_smooth_broken_powerlaw.html)
"""
struct SmoothBrokenPowerlaw{T} <: SpectralModel{T}
    A::T
    γ1::T
    γ2::T
    Eb::T
    m::T
end

"""
    log_sbpl(E, A, γ1, γ2, Eb, m) -> Float64

Log-space evaluation of smooth broken power-law for numerical stability.
"""
function log_sbpl(E, A, γ1, γ2, Eb, m)
    x = E / Eb
    return nm.log(A) - γ1 * nm.log(E) + ((γ1 - γ2) / m) * nm.log(1 + x^m)
end

(m::SmoothBrokenPowerlaw)(E) = exp(log_sbpl(E, m...))

"""
    SmoothBrokenPowerlawFixed{T} <: SpectralModel{T}

Smooth broken power-law with fixed smoothness parameter m=1.

Simplified version with one fewer parameter for easier fitting.
"""
struct SmoothBrokenPowerlawFixed{T} <: SpectralModel{T}
    A::T
    γ1::T
    γ2::T
    Eb::T
end

(m::SmoothBrokenPowerlawFixed)(E) = exp(log_sbpl(E, m..., 1))

# ============================================================================
# Composite Models
# ============================================================================

"""
    TwoStepModel{M1,M2,T}

Two-component spectral model with additive combination.

# Model
```math
f(E) = f_1(E) + f_2(E)
```

Typically used for:
- Low-energy component (f₁): thermal/core population
- High-energy component (f₂): suprathermal/tail population

# Fields
- `model1::M1`: First spectral model
- `model2::M2`: Second spectral model
- `Emin::T`: Transition energy (informational, both models evaluated everywhere)

# Examples
```julia
low_e = PowerLawExpCutoff2(logA=10.0, γ=1.5, logE_c=3.0)
high_e = TransformKappaDistribution(u_A=0.0, u_κ=0.5, u_E_c=0.5)
model = TwoStepModel(low_e, high_e, 30.0)

flux = model.(energies)  # Evaluate at all energies
```
"""
struct TwoStepModel{M1, M2, T}
    model1::M1
    model2::M2
    Emin::T
end

# Backward compatibility: allow indexing
Base.getindex(model::TwoStepModel, i::Integer) = if i == 1
    model.model1
elseif i == 2
    model.model2
elseif i == 3
    model.Emin
else
    throw(BoundsError(model, i))
end

(m::TwoStepModel)(E) = m.model1(E) + m.model2(E)
