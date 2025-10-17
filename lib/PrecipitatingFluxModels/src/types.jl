"""
Core type definitions for precipitating flux models.
"""

"""
    FluxModel

Abstract type for all precipitating flux models.
"""
abstract type FluxModel end

"""
    EmpiricalFluxModel <: FluxModel

Empirical flux model based on statistical analysis of satellite observations.

# Fields
- `parameters::DataFrame`: Statistical parameters organized by spatial bins (MLat, MLT, AE)
- `spatial_bins::NamedTuple`: Bin edges for MLat, MLT, and AE indices
- `energy_range::Tuple{Float64,Float64}`: Valid energy range in keV
- `description::String`: Model description and metadata
- `version::String`: Model version identifier

# Structure of parameters DataFrame
Each row contains:
- `mlat_bin`: Magnetic latitude bin center (degrees)
- `mlt_bin`: Magnetic local time bin center (hours)
- `ae_bin`: AE index bin (categorical)
- `n_samples`: Number of observations in bin
- Statistical parameters (median, Q25, Q75):
  - `κ_*`: Kappa distribution parameter
  - `log_A1_*`, `E_c1_*`, `γ_*`: ExpPow model parameters
  - `log_A2_*`, `E_c2_*`: Kappa model parameters
  - `Emin_*`: Transition energy
  - `J1_*`, `J2_*`: Number flux (ExpPow, Kappa)
  - `JE1_*`, `JE2_*`: Energy flux (ExpPow, Kappa)
"""
struct EmpiricalFluxModel <: FluxModel
    parameters::DataFrame
    spatial_bins::NamedTuple{(:mlat, :mlt, :ae), Tuple{Vector{Float64}, Vector{Float64}, Vector{String}}}
    energy_range::Tuple{Float64, Float64}
    description::String
    version::String
end

"""
    Base.show(io::IO, model::EmpiricalFluxModel)

Pretty printing for EmpiricalFluxModel.
"""
function Base.show(io::IO, model::EmpiricalFluxModel)
    println(io, "EmpiricalFluxModel (v$(model.version))")
    println(io, "  Description: $(model.description)")
    println(io, "  Energy range: $(model.energy_range[1]) - $(model.energy_range[2]) keV")
    println(io, "  Spatial coverage:")
    println(io, "    MLat bins: $(length(model.spatial_bins.mlat))")
    println(io, "    MLT bins: $(length(model.spatial_bins.mlt))")
    println(io, "    AE bins: $(length(model.spatial_bins.ae))")
    println(io, "  Total parameter sets: $(nrow(model.parameters))")
end

"""
    FluxParameters

Container for spectral model parameters at a specific location/condition.
"""
struct FluxParameters
    # ExpPow model (low energy)
    log_A1::Float64
    E_c1::Float64
    γ::Float64

    # Kappa model (high energy)
    log_A2::Float64
    E_c2::Float64
    κ::Float64

    # Transition
    Emin::Float64

    # Integrated fluxes
    J1::Float64  # Number flux (ExpPow)
    J2::Float64  # Number flux (Kappa)
    JE1::Float64  # Energy flux (ExpPow)
    JE2::Float64  # Energy flux (Kappa)

    # Metadata
    mlat::Float64
    mlt::Float64
    ae::Float64
    n_samples::Int
end

"""
    to_spectral_model(params::FluxParameters) -> TwoStepModel

Convert FluxParameters to a TwoStepModel for evaluation.
"""
function to_spectral_model(params::FluxParameters)
    model1 = PowerLawExpCutoff2(params.log_A1, params.γ, log10(params.E_c1))
    model2 = TransformKappaDistribution(params.log_A2, params.κ, log10(params.E_c2))
    return TwoStepModel(model1, model2, params.Emin)
end
