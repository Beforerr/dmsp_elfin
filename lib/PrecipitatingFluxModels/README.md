# PrecipitatingFluxModels.jl

A Julia package for modeling precipitating electron flux using empirical statistical models derived from DMSP and ELFIN observations.

## Features and Roadmap

- 🌍 **Spatial Coverage**: Models flux as a function of magnetic latitude (MLat), magnetic local time (MLT), and geomagnetic activity (AE index)
- 📊 **Statistical Uncertainty**: Provides median values plus 25th and 75th percentiles for all parameters
- 🔬 **Two-Component Spectra**: Combines exponential power-law (thermal) and kappa distribution (suprathermal) models
- [ ] Python bindings (PyCall.jl)

## Installation

```julia
using Pkg

Pkg.add(url="https://github.com/beforerr/dmsp_elfin", subdir="lib/PrecipitatingFluxModels")
```

## Usage Examples

```julia
using PrecipitatingFluxModels

# Load the default model
model = load_model()

# Evaluate flux at specific conditions
# Location: 65° MLat, 6 MLT (dawn sector), moderate activity (AE=150 nT)
j_Efunc = model(; mlat=65.0, mlt=6.0, ae=150.0) # returns a function of energy
flux_10keV = j_Efunc(10.0)

# Get full energy spectrum
energies = 10 .^ range(log10(0.03), log10(1000), length=100)  # 0.03-1000 keV
spectrum = j_Efunc.(energies)

# Compute integrated fluxes
Emin = 0.03
Emax = 100.0
J = n_flux(j_Efunc, Emin, Emax)  # cm⁻² s⁻¹ sr⁻¹
JE = e_flux(j_Efunc, Emin, Emax)  # keV cm⁻² s⁻¹ sr⁻¹
```