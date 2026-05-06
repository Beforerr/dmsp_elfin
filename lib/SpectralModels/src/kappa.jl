# ============================================================================
# Kappa Distribution Models
# ============================================================================

"""
    AbstractKappaDistribution{T} <: SpectralModel{T}

Abstract base type for kappa distribution spectral models.

```math
f(E) = A ⋅ E/E_0 ⋅ [1 + E/(κ⋅E_c)]^{-κ-1}
```

Kappa distributions are similar to the Maxwellian distribution for low and central energies but model suprathermal (high energy) particle populations with power-law tails characterized by the kappa parameter κ.
Setting κ approach infinity leads to the Maxwellian distribution.

# References
- Pierrard & Lazar (2010), doi:10.1007/s11207-010-9640-2
- Vasyliunas (1968), doi:10.1029/JA073i009p02839
- [Plasmapy](https://docs.plasmapy.org/en/stable/api/plasmapy.formulary.distribution.kappa_velocity_1D.html)
"""
abstract type AbstractKappaDistribution{T} <: SpectralModel{T} end

(m::AbstractKappaDistribution)(E) = A(m) * E * (1 + E / (κ(m) * E_c(m)))^(-κ(m) - 1)

function log_eval(m::AbstractKappaDistribution, E)
    return logA(m) + log(E) + log1p(E / (κ(m) * E_c(m))) * (-κ(m) - 1)
end

"""
    KappaDistribution{T} <: AbstractKappaDistribution{T}

Kappa distribution spectral model for suprathermal particle populations.

# Fields
- `A`: Normalization amplitude
- `κ`: Kappa parameter controlling tail behavior (κ > 0)
- `E_c`: Characteristic energy scale (keV)

# Examples
```julia
# Typical magnetospheric suprathermal population
model = KappaDistribution(A=1e4, κ=5.0, E_c=10.0)
flux = model.(10.0:10.0:100.0)
```
"""
@kwdef struct KappaDistribution{T} <: AbstractKappaDistribution{T}
    A::T
    κ::T
    E_c::T
end

KappaDistribution(m::AbstractKappaDistribution) = KappaDistribution(A(m), κ(m), E_c(m))

κ(m::KappaDistribution) = m.κ

# Analytical integrals for number and energy flux
# https://www.wolframalpha.com/input?i=A+x+%281%2Bx%2F%28c+%CE%BA%29%29%5E%28-%CE%BA-1%29
# when integrate from 0 to Inf, the result is (A E_c^2 κ)/(-1 + κ)
function _integral0(m::KappaDistribution, E)
    num = m.A * E_c(m) * m.κ * (E_c(m) + E) * (1 + E / (m.κ * E_c(m)))^(-m.κ)
    den = m.κ - 1
    return -num / den
end

function _integral1(m::KappaDistribution, E)
    κ = m.κ
    c = E_c(m)
    num = A(m) * c * κ * ((κ - 1) * E^2 + 2 * c * E * κ + 2 * c^2 * κ) * (1 + E / (κ * c))^(-κ)
    den = (κ - 1) * (κ - 2)
    return -num / den
end

"""
    TransformKappaDistribution{T} <: AbstractKappaDistribution{T}

Kappa distribution with bounded parameter transformations for optimization.

Uses TransformVariables.jl to map unconstrained optimization parameters to
valid physical ranges:
- `A`: ℝ → ℝ⁺ (positive)
- `κ`: ℝ → [1, 20] (physically meaningful range)
- `E_c`: ℝ → [0.01, 100] keV (typical particle energies)

# Fields
- `u_A`: Transformed amplitude parameter
- `u_κ`: Transformed kappa parameter
- `u_E_c`: Transformed energy parameter

# Examples
```julia
# Use with unconstrained optimization
init_params = [0.0, 0.5, 0.5]  # Transformed space
model = TransformKappaDistribution(u_A=init_params[1], u_κ=init_params[2], u_E_c=init_params[3])

# Access physical parameters
A_phys = A(model)      # Gets transformed A > 0
κ_phys = κ(model)      # Gets transformed κ ∈ [1, 20]
E_c_phys = E_c(model)  # Gets transformed E_c ∈ [0.01, 100]
```
"""
@kwdef struct TransformKappaDistribution{T} <: AbstractKappaDistribution{T}
    u_A::T
    u_κ::T
    u_E_c::T
end

# Transform specifications
const as_A = asℝ₊                    # A > 0
const as_κ = as(Real, 1, 20)         # κ ∈ [1, 20]
const as_Ec = as(Real, 0.01, 100)    # E_c ∈ [0.01, 100] keV

@inline A(m::TransformKappaDistribution) = transform(as_A, m.u_A)
@inline κ(m::TransformKappaDistribution) = transform(as_κ, m.u_κ)
@inline E_c(m::TransformKappaDistribution) = transform(as_Ec, m.u_E_c)
logA(m::TransformKappaDistribution) = log(transform(as_A, m.u_A))

TransformKappaDistribution(m::AbstractKappaDistribution) = begin
    u1 = inverse(as_A, A(m))
    u2 = inverse(as_κ, κ(m))
    u3 = inverse(as_Ec, min(E_c(m), 100))
    return TransformKappaDistribution(u1, u2, u3)
end

"""
    KappaDistribution2{T} <: AbstractKappaDistribution{T}

Kappa distribution with log-transformed parameters `logA`, `logκ`, and `logE_c`.

Where `A = exp(logA)`, `κ = exp(logκ)`, and `E_c = exp(logE_c)`.

This parametrization ensures all physical parameters remain positive during optimization.

# Examples
```julia
# Create model with log-space parameters
model = KappaDistribution2(logA=9.2, logκ=1.6, logE_c=2.3)

# Evaluate
flux = model(50.0)

# Access physical parameters
A_phys = A(model)      # exp(9.2)
κ_phys = κ(model)      # exp(1.6) ≈ 5.0
E_c_phys = E_c(model)  # exp(2.3) ≈ 10.0 keV
```
"""
struct KappaDistribution2{T} <: AbstractKappaDistribution{T}
    logA::T
    logκ::T
    logE_c::T
end

κ(m::KappaDistribution2) = exp(m.logκ)

KappaDistribution2(m::AbstractKappaDistribution) = KappaDistribution2(logA(m), log(κ(m)), logE_c(m))


_integral0(m::AbstractKappaDistribution, E) = _integral0(KappaDistribution(m), E)
_integral1(m::AbstractKappaDistribution, E) = _integral1(KappaDistribution(m), E)

# ============================================================================
# Display and Utilities
# ============================================================================

"""
    math_show(model::AbstractKappaDistribution; sigdigits=2) -> LaTeXString

Return LaTeX representation of kappa distribution with parameter values.

# Examples
```julia
model = KappaDistribution(A=1e4, κ=5.0, E_c=10.0)
math_show(model)  # L"A = 1.0e4, κ = 5.0, E_c = 10"
```
"""
function math_show(m::AbstractKappaDistribution; sigdigits = 2)
    _A = round(A(m); sigdigits)
    _κ = round(κ(m); sigdigits)
    _E_c = round(E_c(m); sigdigits)
    return L"""
    $A_κ ⋅ E/E_0 ⋅ (1 + E/κ E_κ)^{-κ-1}$
    \\
    $A_κ$:%$_A, $κ$:%$_κ, $E_κ$:%$_E_c
    """
end
