"""
    SpectralModel{T}

Abstract base type for all spectral models.

Features:
- Callable: `model(E)` evaluates flux at energy E
- Parameter extraction: `raw_vec(model)`
"""
abstract type SpectralModel{T} end

function Base.iterate(m::T, state = 1) where {T <: SpectralModel}
    return state > fieldcount(T) ? nothing : (getfield(m, state), state + 1)
end

Base.broadcastable(o::SpectralModel) = Ref(o)

## Accessors
A(m) = hasproperty(m, :A) ? m.A : exp(m.logA)
logA(m) = hasproperty(m, :logA) ? m.logA : log(A(m))
E_c(m) = hasproperty(m, :E_c) ? m.E_c : exp(m.logE_c)
logE_c(m) = hasproperty(m, :logE_c) ? m.logE_c : log(E_c(m))

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
    Tbare = Base.typename(T).wrapper  # strip type param so Dual numbers are accepted
    return quote
        @assert length(u) == $n "Expected $($n) parameters, got $(length(u))"
        @inbounds $(Tbare)($(args...))
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
