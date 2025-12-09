"""
Core type definitions for precipitating flux models.
"""

abstract type FluxModel end

"""
    EmpiricalFluxModel <: FluxModel

Empirical flux model based on statistical analysis of satellite observations.

# Fields
- `parameters::DataFrame`: Statistical parameters organized by spatial bins (MLat, MLT, AE)
- `spatial_bins::NamedTuple`: Bin edges for MLat, MLT, and AE indices
- `energy_range::Tuple{Float64,Float64}`: Valid energy range in keV
- `description`: Model description and metadata
- `version`: Model version identifier
"""
struct EmpiricalFluxModel{D, S} <: FluxModel
    parameters::D
    spatial_bins::S
    energy_range::Tuple{Float64, Float64}
    description::String
    version::String
end

function Base.show(io::IO, model::EmpiricalFluxModel)
    println(io, "EmpiricalFluxModel (v$(model.version))")
    println(io, "  Description: $(model.description)")
    println(io, "  Energy range: $(model.energy_range[1]) - $(model.energy_range[2]) keV")
    println(io, "  Spatial coverage:")
    println(io, "    MLat bins: $(length(model.spatial_bins.mlat))")
    println(io, "    MLT bins: $(length(model.spatial_bins.mlt))")
    println(io, "    AE bins: $(length(model.spatial_bins.ae))")
    return println(io, "  Total models: $(nrow(model.parameters))")
end

function (m::EmpiricalFluxModel)(; mlat = 65.0, mlt = 6.0, ae = 150.0)
    return get_model(m; mlat, mlt, ae)
end

function (m::EmpiricalFluxModel)(E)
    return mean(x -> x(E), m.parameters.model)
end
