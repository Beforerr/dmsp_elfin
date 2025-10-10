# ============================================================================
# Kappa Distribution Models
# ============================================================================

"""
    AbstractKappaDistribution{T} <: SpectralModel{T}

Abstract base type for kappa distribution spectral models.

```math
f(E) = A ⋅ E ⋅ [1 + E/(κ⋅E_c)]^{-κ-1}
```

Kappa distributions are similar to the Maxwellian distribution for low and central energies but model suprathermal (high energy) particle populations with power-law tails characterized by the kappa parameter κ.
Setting κ approach infinity leads to the Maxwellian distribution.

# References
- Pierrard & Lazar (2010), doi:10.1007/s11207-010-9640-2
- Vasyliunas (1968), doi:10.1029/JA073i009p02839
- [Plasmapy](https://docs.plasmapy.org/en/stable/api/plasmapy.formulary.distribution.kappa_velocity_1D.html)
"""
abstract type AbstractKappaDistribution{T} <: SpectralModel{T} end

"""
    KappaDistribution{T} <: AbstractKappaDistribution{T}

Kappa distribution spectral model for suprathermal particle populations.

# Fields
- `A::T`: Normalization amplitude
- `κ::T`: Kappa parameter controlling tail behavior (κ > 0)
- `E_c::T`: Characteristic energy scale (keV)

# Examples
```julia
# Typical magnetospheric suprathermal population
model = KappaDistribution(A=1e4, κ=5.0, E_c=10.0)
flux = model.(10.0:10.0:100.0)

# Compute integrated quantities
n = n_flux(model, 10.0, 100.0)  # Total particle count
e = e_flux(model, 10.0, 100.0)  # Total energy
```
"""
@kwdef struct KappaDistribution{T} <: AbstractKappaDistribution{T}
    A::T
    κ::T
    E_c::T
end

(m::KappaDistribution)(E) = m.A * E * (1 + E / (m.κ * m.E_c))^(-m.κ - 1)

@inline function log_eval(m::KappaDistribution, E)
    return nm.log(m.A) + log(E) + nm.log(1 + E / (m.κ * m.E_c)) * (-m.κ - 1)
end

# Accessor functions
A(m::KappaDistribution) = m.A
E_c(m::KappaDistribution) = m.E_c
κ(m::KappaDistribution) = m.κ
logA(m::KappaDistribution) = log(m.A)

# Analytical integrals for number and energy flux
# https://www.wolframalpha.com/input?i=A+x+%281%2Bx%2F%28c+%CE%BA%29%29%5E%28-%CE%BA-1%29
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

# ============================================================================
# Transformed Kappa Distributions
# ============================================================================

"""
    TransformKappaDistribution{T} <: AbstractKappaDistribution{T}

Kappa distribution with bounded parameter transformations for optimization.

Uses TransformVariables.jl to map unconstrained optimization parameters to
valid physical ranges:
- `A`: ℝ → ℝ⁺ (positive)
- `κ`: ℝ → [1, 20] (physically meaningful range)
- `E_c`: ℝ → [0.01, 100] keV (typical particle energies)

# Fields
- `u_A::T`: Transformed amplitude parameter
- `u_κ::T`: Transformed kappa parameter
- `u_E_c::T`: Transformed energy parameter

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

# Accessor functions with transformations
@inline A(m::TransformKappaDistribution) = transform(as_A, m.u_A)
@inline κ(m::TransformKappaDistribution) = transform(as_κ, m.u_κ)
@inline E_c(m::TransformKappaDistribution) = transform(as_Ec, m.u_E_c)

# Evaluation
(m::TransformKappaDistribution)(E) = KappaDistribution(m)(E)
KappaDistribution(m::TransformKappaDistribution) = KappaDistribution(A(m), κ(m), E_c(m))

@inline function log_eval(m::TransformKappaDistribution, E)
    return log_eval(KappaDistribution(m), E)
end

# ============================================================================
# Log-Transformed Kappa Distribution
# ============================================================================

"""
    KappaDistribution2{T} <: AbstractKappaDistribution{T}

Kappa distribution with log-transformed parameters for unconstrained optimization.

# Model
```math
f(E) = A ⋅ E ⋅ [1 + E/(κ⋅E_c)]^{-κ-1}
```

where `A = exp(logA)`, `κ = exp(logκ)`, and `E_c = exp(logE_c)`.

# Fields
- `logA::T`: Log of normalization amplitude
- `logκ::T`: Log of kappa parameter
- `logE_c::T`: Log of characteristic energy

This parametrization ensures all physical parameters remain positive during
unconstrained optimization over ℝ³.

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

(m::KappaDistribution2)(E) = A(m) * E * (1 + E / (κ(m) * E_c(m)))^(-κ(m) - 1)

function log_eval(m::KappaDistribution2, E)
    return m.logA + log(E) + log(1 + E / (κ(m) * E_c(m))) * (-κ(m) - 1)
end

# Accessor functions
A(m::KappaDistribution2) = exp(m.logA)
logA(m::KappaDistribution2) = m.logA
κ(m::KappaDistribution2) = exp(m.logκ)
E_c(m::KappaDistribution2) = exp(m.logE_c)

# Convert to base KappaDistribution for integrals
KappaDistribution(m::KappaDistribution2) = KappaDistribution(A(m), κ(m), E_c(m))
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
    $A ⋅ E ⋅ \left(1 + E/κ E_c\right)^{-κ-1}$
    \\
    A:%$_A, κ:%$_κ, E_c:%$_E_c
    """
end

"""
    log_jacobian(model::KappaDistribution, E) -> Vector{Float64}

Compute the Jacobian of log(f(E)) with respect to model parameters [A, κ, E_c].

Used for gradient-based optimization and uncertainty quantification.

# Returns
Vector [∂log(f)/∂A, ∂log(f)/∂κ, ∂log(f)/∂E_c]
"""
function log_jacobian(m::KappaDistribution, E)
    A_val = m.A
    κ_val = m.κ
    E_c_val = m.E_c

    u = E / (κ_val * E_c_val)  # Normalized energy
    g = nm.log(1 + u)

    # Partial derivatives
    dA = 1 / A_val
    dκ = -g + (κ_val + 1) * E / (κ_val * (κ_val * E_c_val + E))
    dEc = (κ_val + 1) * E / (E_c_val * (κ_val * E_c_val + E))

    return [dA, dκ, dEc]
end
