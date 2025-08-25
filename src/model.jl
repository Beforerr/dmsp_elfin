using NaNMath
import NaNMath as nm

# https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/plot_exp_cutoff_powerlaw.html

"""
f(E) = A * E^(-γ) * exp(-E/E_c)
"""
struct PowerLawExp{T}
    A::T
    γ::T
    E_c::T
end

(f::PowerLawExp)(E) = f.A * E^(-f.γ) * exp(-E / f.E_c)

# https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/plot_smooth_broken_powerlaw.html
struct SmoothBrokenPowerlaw{T}
    A::T
    γ1::T
    γ2::T
    Eb::T
    m::T
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

(f::SmoothBrokenPowerlaw)(E) = sbpl(E, f.A, f.γ1, f.γ2, f.Eb, f.m)