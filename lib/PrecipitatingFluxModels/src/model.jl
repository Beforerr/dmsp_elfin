"""
    get_model(model::EmpiricalFluxModel; mlat=nothing, mlt=nothing, ae=nothing)

Query subset of models matching spatial bin criteria.

Returns an `EmpiricalFluxModel` containing only rows within the specified bins.
The returned model can be called with energy `E` to compute mean flux.

# Example
```julia
model = load_model()
subset = get_model(model; mlat=70, ae=200)  # Models in mlat bin containing 70°, AE bin containing 200
flux = subset(50.0)  # Mean flux at 50 keV across matching models
```
"""
function get_model(m; mlat = nothing, mlt = nothing, ae = nothing)
    df = m.parameters
    bins = m.spatial_bins

    if !isnothing(mlat)
        lo, hi = _find_bin_edges(abs(mlat), bins.mlat)
        df = @subset(df, lo .<= abs.(:mlat) .< hi; view = true)
    end
    if !isnothing(mlt)
        lo, hi = _find_bin_edges(mlt, bins.mlt)
        df = @subset(df, lo .<= :mlt .< hi; view = true)
    end
    if !isnothing(ae)
        lo, hi = _find_bin_edges(ae, bins.ae)
        df = @subset(df, lo .<= :ae .< hi; view = true)
    end

    return EmpiricalFluxModel(df, bins, m.energy_range, m.description, m.version)
end

# Find bin edges [lo, hi) containing val
function _find_bin_edges(val, bins)
    for i in 1:(length(bins) - 1)
        bins[i] <= val < bins[i + 1] && return (bins[i], bins[i + 1])
    end
    # Handle edge case: val >= last bin edge
    return (bins[end - 1], bins[end])
end

"""
    n_flux(model::EmpiricalFluxModel, Emin, Emax)

Compute mean number flux across all models in the collection.
"""
function n_flux(m::EmpiricalFluxModel, Emin, Emax)
    return mean(x -> SpectralModels.n_flux(x, Emin, Emax), m.parameters.model)
end

"""
    e_flux(model::EmpiricalFluxModel, Emin, Emax)

Compute mean energy flux across all models in the collection.
"""
function e_flux(m::EmpiricalFluxModel, Emin, Emax)
    return mean(x -> SpectralModels.e_flux(x, Emin, Emax), m.parameters.model)
end
