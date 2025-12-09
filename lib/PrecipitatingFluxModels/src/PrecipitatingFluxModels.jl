module PrecipitatingFluxModels

using Statistics
using DataFrames, DataFramesMeta
using SpectralModels
using JLD2

# Re-export from SpectralModels for convenience
export TwoStepModel, PowerLawExpCutoff2, TransformKappaDistribution
export n_flux, e_flux

# Package exports
export EmpiricalFluxModel
export load_model, get_model

include("types.jl")
include("io.jl")
include("model.jl")

end
