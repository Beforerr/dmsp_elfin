# ============================================================================
# Power-Law Models
# ============================================================================

"""
    PowerLaw{T} <: SpectralModel{T}

Simple power-law spectral model.

# Model
```math
f(E) = A ⋅ E^{-γ}
```

# Fields
- `A::T`: Normalization amplitude
- `γ::T`: Spectral index (positive for falling spectra)

# Examples
```julia
model = PowerLaw(1.5e5, 2.1)
flux = model(50.0)  # Evaluate at 50 keV
```
"""
struct PowerLaw{T} <: SpectralModel{T}
    A::T
    γ::T
end

(f::PowerLaw)(E) = f.A * E^(-f.γ)

"""
    PowerLawExpCutoff{T} <: SpectralModel{T}

Power-law with exponential cutoff.

# Model
```math
f(E) = A ⋅ E^{-γ} ⋅ \\exp(-E/E_c)
```

# Fields
- `A::T`: Normalization amplitude
- `γ::T`: Spectral index
- `E_c::T`: Cutoff energy (must be positive)

# Examples
```julia
model = PowerLawExpCutoff(A=1e5, γ=2.0, E_c=50.0)
flux = model.(10.0:10.0:100.0)
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
    PowerLawExpCutoff2{T} <: SpectralModel{T}

Power-law with exponential cutoff using log-transformed parameters.

# Model
```math
f(E) = A ⋅ E^{-γ} ⋅ \\exp(-E/E_c)
```

where `A = exp(logA)` and `E_c = exp(logE_c)`.

# Fields
- `logA::T`: Log of normalization amplitude
- `γ::T`: Spectral index
- `logE_c::T`: Log of cutoff energy

This parametrization ensures A > 0 and E_c > 0 during unconstrained optimization.

# Examples
```julia
model = PowerLawExpCutoff2(logA=11.5, γ=2.0, logE_c=3.9)
flux = model(50.0)
```
"""
@kwdef struct PowerLawExpCutoff2{T} <: SpectralModel{T}
    logA::T
    γ::T
    logE_c::T
end

(m::PowerLawExpCutoff2)(E) = A(m) * E^(-m.γ) * exp(-E / E_c(m))

A(m::PowerLawExpCutoff2) = exp(m.logA)
logA(m::PowerLawExpCutoff2) = m.logA
E_c(m::PowerLawExpCutoff2) = exp(m.logE_c)
log_eval(m::PowerLawExpCutoff2, E) = m.logA - m.γ * log(E) - E / E_c(m)

_clean_text(x) = replace(x, "+0" => "", "e+" => "e", "e0" => "e")

function math_show(m::PowerLawExpCutoff2; sigdigits = 2)
    # _A = @sprintf("%.2g", exp(m.logA))
    # _E_c = @sprintf("%.2g", exp(m.logE_c))
    # _γ = @sprintf("%.2g", m.γ)

    _A = round(A(m); sigdigits)
    _E_c = @sprintf("%.2g", E_c(m)) |> _clean_text
    _γ = @sprintf("%.2g", m.γ) |> _clean_text

    return L"""
    $A ⋅ E^{-γ} ⋅ e^{-E/E_c}$
    \\
    A:%$_A, γ:%$_γ, E_c:%$_E_c
    """
end

# Analytical integrals for PowerLawExpCutoff2
# https://www.wolframalpha.com/input?i=A+x%5E%28-+%CE%B3%29+exp%28-x%2Fc%29
_integral0(m::PowerLawExpCutoff2, E) = -A(m) * E_c(m)^(1 - m.γ) * gamma(1 - m.γ, E / E_c(m))
_integral1(m::PowerLawExpCutoff2, E) = -A(m) * E_c(m)^(2 - m.γ) * gamma(2 - m.γ, E / E_c(m))
