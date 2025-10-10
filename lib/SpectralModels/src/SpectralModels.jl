"""
    SpectralModels

Spectral model fitting for particle flux data.

This module provides spectral models (power-law, exponential cutoff, kappa distributions)
and fitting routines for particle flux measurements.

# Main Types

## Spectral Models
- `PowerLaw`: Simple power-law f(E) = A·E^(-γ)
- `PowerLawExpCutoff`: Power-law with exponential cutoff f(E) = A·E^(-γ)·exp(-E/E_c)
- `PowerLawExpCutoff2`: Log-transformed version for unconstrained optimization
- `KappaDistribution`: Kappa distribution f(E) = A·E·(1 + E/(κ·E_c))^(-κ-1)
- `KappaDistribution2`: Log-transformed kappa distribution
- `TransformKappaDistribution`: Bounded kappa distribution for optimization
- `SmoothBrokenPowerlaw`: Smooth broken power-law model
- `SmoothBrokenPowerlawFixed`: Fixed smoothness parameter version
- `TwoStepModel`: Combined two-component model

# Main Functions

## Fitting
- `fit(ModelType, energies, flux)`: Fit spectral model to data
- `init_guess(ModelType, energies, flux)`: Generate initial parameter guess
- `fit_flux_two_step`: Fit two-component model with transition energy
- `fit_two_flux`: Fit combined flux from two instruments

## Flux Integration
- `n_flux(model, Emin, Emax)`: Compute number flux in energy range
- `e_flux(model, Emin, Emax)`: Compute energy flux in energy range

## Utilities
- `paramcount(model)`: Number of model parameters
- `raw_vec(model)`: Extract parameter vector
- `log_eval(model, E)`: Log-space evaluation
- `math_show(model)`: LaTeX representation of model

# Examples

```julia
using SpectralModels

# Simple power-law fit
energies = [10.0, 20.0, 50.0, 100.0]  # keV
flux = [1e4, 5e3, 1e3, 200.0]  # particles/(cm²·s·sr·keV)

model = fit(PowerLaw, energies, flux)
# PowerLaw(A=1.5e5, γ=2.1)

# Evaluate model at new energies
new_energies = 15.0:5.0:95.0
predicted_flux = model.(new_energies)

# Power-law with exponential cutoff
model2 = fit(PowerLawExpCutoff, energies, flux)

# Compute integrated flux
total_flux = n_flux(model2, 10.0, 100.0)
energy_flux = e_flux(model2, 10.0, 100.0)

# Kappa distribution for suprathermal populations
kappa_model = fit(KappaDistribution, energies, flux)

# Two-component model
two_step = fit(TwoStepModel{PowerLawExpCutoff2, TransformKappaDistribution},
               flux, energies)
```

# References

- Gammapy spectral models: https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/
- SpectralFitting.jl: https://fjebaker.github.io/SpectralFitting.jl/dev/
- sherpa: https://github.com/sherpa/sherpa : a Python package for modeling and fitting data.
"""
module SpectralModels

# Standard library
using Printf
using Statistics: mean

# External dependencies
using LaTeXStrings
using TransformVariables: transform, as, inverse, asℝ₊
using NaNMath
import NaNMath as nm
using LsqFit
using StaticArrays
using NonlinearSolve
using CurveFit
using SpecialFunctions: gamma
import Base: iterate, show

export SpectralModel, AbstractKappaDistribution
export PowerLaw, PowerLawExpCutoff, PowerLawExpCutoff2
export KappaDistribution, KappaDistribution2, TransformKappaDistribution
export SmoothBrokenPowerlaw, SmoothBrokenPowerlawFixed
export TwoStepModel

export fit, init_guess
export fit_flux_two_step, fit_two_flux
export sciml_log_fit, LsqFit_log_fit

export paramcount, raw_vec, log_eval
export n_flux, e_flux
export A, logA, E_c, κ, math_show

include("model.jl")
include("kappa.jl")
include("powerlaw.jl")
include("fit.jl")
include("sciml.jl")

end
