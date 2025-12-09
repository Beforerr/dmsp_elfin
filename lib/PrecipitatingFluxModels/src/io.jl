"""
Input/output functions for model serialization.
"""

const DEFAULT_PATH = joinpath(pkgdir(PrecipitatingFluxModels), "release", "default_model.jld2")

save(model, fname = DEFAULT_PATH) = @save fname model

function load_model(fname = DEFAULT_PATH)
    @load fname model
    return model
end

"""
    create_model_from_stats(df;
        mlat_bins, mlt_bins, ae_bins, energy_range, description, version
    )

Create an EmpiricalFluxModel from statistical analysis results `df`.

# Example
```julia
model = create_model_from_stats(df;
    mlat_bins=50:1:80,
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
