# ============================================================================
# Power-Law Models
# ============================================================================
abstract type AbstractPowerLaw{T} <: SpectralModel{T} end

const PowerLawExpCutoffMath = """
```math
f(E) = A ⋅ (E/E_0)^{-γ} ⋅ \\exp(-E/E_c)
```
"""

abstract type AbstractPowerLawExpCutoff{T} <: AbstractPowerLaw{T} end

(m::AbstractPowerLaw)(E) = A(m) * E^(-m.γ)
(m::AbstractPowerLawExpCutoff)(E) = A(m) * E^(-m.γ) * exp(-E / E_c(m))

log_eval(m::AbstractPowerLawExpCutoff, E) = logA(m) - m.γ * log(E) - E / E_c(m)

"""
    PowerLaw{T} <: SpectralModel{T}

Simple power-law spectral model.

# Model
```math
f(E) = A ⋅ E^{-γ}
```

# Fields
- `A`: Normalization amplitude
- `γ`: Spectral index (positive for falling spectra)

# Examples
```julia
model = PowerLaw(1.5e5, 2.1)
flux = model(50.0)  # Evaluate at 50 keV
```
"""
struct PowerLaw{T} <: AbstractPowerLaw{T}
    A::T
    γ::T
end

"""
    PowerLawExpCutoff{T} <: AbstractPowerLaw{T}

Power-law with exponential cutoff.

# Fields
- `A`: Normalization amplitude
- `γ`: Spectral index
- `E_c`: Cutoff energy (must be positive)

# Examples
```julia
model = PowerLawExpCutoff(A=1e5, γ=2.0, E_c=50.0)
flux = model.(10.0:10.0:100.0)
```
"""
@kwdef struct PowerLawExpCutoff{T} <: AbstractPowerLawExpCutoff{T}
    A::T
    γ::T
    E_c::T
end

"""
    PowerLawExpCutoff2{T}

Power-law with exponential cutoff using log-transformed parameters `logA` and `logE_c`.

$PowerLawExpCutoffMath

where `A = exp(logA)` and `E_c = exp(logE_c)`.

This parametrization ensures A > 0 and E_c > 0 during optimization.

# Examples
```julia
model = PowerLawExpCutoff2(logA=11.5, γ=2.0, logE_c=3.9)
flux = model(50.0)
```
"""
@kwdef struct PowerLawExpCutoff2{T} <: AbstractPowerLawExpCutoff{T}
    logA::T
    γ::T
    logE_c::T
end

PowerLawExpCutoff2(m::AbstractPowerLawExpCutoff) = PowerLawExpCutoff2(logA(m), m.γ, logE_c(m))

_clean_text(x) = replace(x, "+0" => "", "e+" => "e", "e0" => "e")

function math_show(m::PowerLawExpCutoff2; sigdigits = 2)
    # _A = @sprintf("%.2g", exp(m.logA))
    # _E_c = @sprintf("%.2g", exp(m.logE_c))
    # _γ = @sprintf("%.2g", m.γ)

    _A = round(A(m); sigdigits)
    _E_c = @sprintf("%.2g", E_c(m)) |> _clean_text
    _γ = @sprintf("%.2g", m.γ) |> _clean_text

    return L"""
    $A_{EP} ⋅ (E/E_0)^{-γ} ⋅ \exp(-E/E_{EP})$
    \\
    $A_{EP}$:%$_A, $γ$:%$_γ, $E_{EP}$:%$_E_c
    """
end

# Analytical integrals for PowerLawExpCutoff2
# https://www.wolframalpha.com/input?i=A+x%5E%28-+%CE%B3%29+exp%28-x%2Fc%29
_integral0(m::PowerLawExpCutoff2, E) = -A(m) * E_c(m)^(1 - m.γ) * gamma(1 - m.γ, E / E_c(m))
_integral1(m::PowerLawExpCutoff2, E) = -A(m) * E_c(m)^(2 - m.γ) * gamma(2 - m.γ, E / E_c(m))

function unbound_fit(M::Type{<:PowerLawExpCutoff}, E, y)
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
    return M(A, γ, E_c)
end
