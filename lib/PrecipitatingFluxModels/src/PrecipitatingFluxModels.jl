"""
    PrecipitatingFluxModels

A Julia package for modeling precipitating electron flux using empirical statistical models
derived from DMSP and ELFIN satellite observations.

# Main Features
- Evaluate precipitating flux spectra as a function of MLat, MLT, and AE index
- Support for both individual spectral components and combined models
- Statistical parameters (median, percentiles) for model uncertainty quantification

# Exported Functions
- `evaluate_flux`: Evaluate flux at specific geophysical conditions
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
export load_model
export evaluate_flux
export flux_parameters
export number_flux, energy_flux

include("types.jl")
include("model.jl")
include("evaluation.jl")
include("io.jl")

end
