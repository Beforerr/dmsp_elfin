"""
    SmoothBrokenPowerlaw{T} <: SpectralModel{T}

Smooth broken power-law model with variable smoothness.

# Model
```math
f(E) = A ⋅ E^{-γ_1} ⋅ [1 + (E/E_b)^m]^{(γ_2 - γ_1)/m}
```

# Fields
- `A`: Normalization amplitude
- `γ1`: Low-energy spectral index
- `γ2`: High-energy spectral index
- `Eb`: Break energy
- `m`: Smoothness parameter (m=1 gives smooth transition)

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
    return log(A) - γ1 * log(E) + ((γ1 - γ2) / m) * log(1 + x^m)
end

(m::SmoothBrokenPowerlaw)(E) = exp(log_sbpl(E, m...))
