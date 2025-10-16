"""
Input/output functions for model serialization.
"""

save(fname, model) = @save fname model

function load_model(fname)
    @load fname model
    return model
end

function load_model()
    pkg_dir = pkgdir(PrecipitatingFluxModels)
    default_model = joinpath(pkg_dir, "release", "default_model.jld2")
    return load_model(default_model)
end

"""
    create_model_from_stats(df; mlat_bins, mlt_bins, ae_bins,
        energy_range=(0.03, 1000.0), description="Empirical flux model", version="1.0.0"
    )

Create an EmpiricalFluxModel from statistical analysis results `df`.

# Arguments
- `mlat_bins`: Vector of MLat bin edges
- `mlt_bins`: Vector of MLT bin edges
- `ae_bins`: Vector of AE bin edges
- `energy_range`: Valid energy range in keV

# Example
```julia
model = create_model_from_stats(df;
    mlat_bins=50:1:80,
    mlt_bins=0:2:22,
    ae_bins=[0, 100, 300, Inf],
)
```
"""
function create_model_from_stats(
        df;
        mlat_bins = 52.5:1:80.5,
        mlt_bins = 0:1:24,
        ae_bins = [0, 100, 300, Inf],
        energy_range = (0.03, 1000.0),
        description = "DMSP-ELFIN 2020-2022 model",
        version = "1.0.0"
    )

    # Clean the dataframe
    sdf = dropmissing!(@select(df, :model, :mlat, :mlt = :mlt_elx, :ae = :maxAE))

    # Create spatial bins named tuple
    spatial_bins = (
        mlat = collect(Float64, mlat_bins),
        mlt = collect(Float64, mlt_bins),
        ae = collect(Float64, ae_bins),
    )

    return EmpiricalFluxModel(
        sdf,
        spatial_bins,
        Tuple(Float64.(energy_range)),
        description,
        version
    )
end
