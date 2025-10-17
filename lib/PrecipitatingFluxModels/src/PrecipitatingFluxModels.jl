"""
    PrecipitatingFluxModels

A Julia package for modeling precipitating electron flux using empirical statistical models
derived from DMSP and ELFIN satellite observations.

# Main Features
- Evaluate precipitating flux spectra as a function of MLat, MLT, and AE index
- Support for both individual spectral components and combined models
- Statistical parameters (median, percentiles) for model uncertainty quantification

# Quick Start
```julia
using PrecipitatingFluxModels

# Load the default model
model = load_flux_model()

# Evaluate flux at specific conditions
flux = evaluate_flux(model, mlat=65.0, mlt=6.0, ae=150.0)

# Get energy spectrum
energies = 10 .^ range(log10(0.03), log10(1000), length=100)  # keV
spectrum = flux_spectrum(model, energies; mlat=65.0, mlt=6.0, ae=150.0)
```

# Exported Functions
- `load_flux_model`: Load empirical flux model from data
- `evaluate_flux`: Evaluate flux at specific geophysical conditions
- `flux_spectrum`: Compute full energy spectrum
- `flux_parameters`: Get model parameters for given conditions
- `number_flux`: Compute integrated number flux
- `energy_flux`: Compute integrated energy flux
"""
module PrecipitatingFluxModels

using Statistics
using DataFrames, DataFramesMeta
using SpectralModels
using JLD2

# Re-export from SpectralModels for convenience
export TwoStepModel, PowerLawExpCutoff2, TransformKappaDistribution
export n_flux, e_flux

# Package exports
export FluxModel, EmpiricalFluxModel
export load_flux_model, save_flux_model
export evaluate_flux, flux_spectrum
export flux_parameters
export number_flux, energy_flux
export get_ae_bin, interpolate_parameters

include("types.jl")
include("model.jl")
include("evaluation.jl")
include("io.jl")

end
