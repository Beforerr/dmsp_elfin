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

e_flux(m::TwoStepModel, Emin, Emax) =
    e_flux(m.model1, Emin, Emax) +
    e_flux(m.model2, Emin, Emax)

n_flux(m::TwoStepModel, Emin, Emax) =
    n_flux(m.model1, Emin, Emax) +
    n_flux(m.model2, Emin, Emax)
