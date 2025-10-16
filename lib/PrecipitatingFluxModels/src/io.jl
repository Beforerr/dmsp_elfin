"""
Input/output functions for model serialization.
"""

"""
    save_flux_model(filename::String, model::EmpiricalFluxModel)

Save empirical flux model to a JLD2 file.

# Example
```julia
save_flux_model("flux_model_v1.jld2", model)
```
"""
function save_flux_model(filename::String, model::EmpiricalFluxModel)
    jldsave(filename;
        parameters=model.parameters,
        spatial_bins=model.spatial_bins,
        energy_range=model.energy_range,
        description=model.description,
        version=model.version
    )
    @info "Model saved to $filename"
end

"""
    load_flux_model(filename::String) -> EmpiricalFluxModel

Load empirical flux model from a JLD2 file.

# Example
```julia
model = load_flux_model("flux_model_v1.jld2")
```
"""
function load_flux_model(filename::String)
    data = load(filename)
    return EmpiricalFluxModel(
        data["parameters"],
        data["spatial_bins"],
        data["energy_range"],
        data["description"],
        data["version"]
    )
end

"""
    load_flux_model() -> EmpiricalFluxModel

Load the default empirical flux model included with the package.

# Example
```julia
model = load_flux_model()
```
"""
function load_flux_model()
    # Look for default model in package data directory
    pkg_dir = pkgdir(PrecipitatingFluxModels)
    default_model = joinpath(pkg_dir, "data", "default_model.jld2")

    if !isfile(default_model)
        error("""
        Default model not found at: $default_model

        To create a default model:
        1. Prepare your statistical parameters DataFrame
        2. Create an EmpiricalFluxModel
        3. Save it using save_flux_model()
        """)
    end

    return load_flux_model(default_model)
end

"""
    create_model_from_stats(df::DataFrame;
                           mlat_bins, mlt_bins, ae_bins,
                           energy_range=(0.03, 1000.0),
                           description="Empirical flux model",
                           version="1.0.0") -> EmpiricalFluxModel

Create an EmpiricalFluxModel from statistical analysis results.

# Arguments
- `df`: DataFrame with columns for each parameter (with _median, _q25, _q75 suffixes)
- `mlat_bins`: Vector of MLat bin centers
- `mlt_bins`: Vector of MLT bin centers
- `ae_bins`: Vector of AE bin labels (e.g., ["[0, 100)", "[100, 300)"])
- `energy_range`: Valid energy range in keV
- `description`: Model description
- `version`: Model version string

# Required DataFrame Columns
The DataFrame must contain columns for spatial binning and statistical parameters:
- `mlat_bin`: Magnetic latitude bin center
- `mlt_bin`: Magnetic local time bin center
- `ae_bin`: AE index bin (categorical)
- `n_samples`: Number of observations
- For each parameter X: `X_median`, `X_q25`, `X_q75`
  where X ∈ {κ, log_A1, E_c1, γ, log_A2, E_c2, Emin, J1, J2, JE1, JE2}

# Example
```julia
# Prepare statistics from your analysis
df = compute_spatial_statistics(raw_data)

# Create model
model = create_model_from_stats(df;
    mlat_bins=50:1:80,
    mlt_bins=0:2:22,
    ae_bins=["[0, 100)", "[100, 300)", "[300, Inf)"],
    description="DMSP-ELFIN 2020-2022 model",
    version="1.0.0"
)

# Save for distribution
save_flux_model("dmsp_elfin_model_v1.jld2", model)
```
"""
function create_model_from_stats(df::DataFrame;
                                mlat_bins::AbstractVector,
                                mlt_bins::AbstractVector,
                                ae_bins::AbstractVector{String},
                                energy_range::Tuple{Real,Real}=(0.03, 1000.0),
                                description::String="Empirical flux model",
                                version::String="1.0.0")

    # Validate required columns
    required_base_cols = [:mlat_bin, :mlt_bin, :ae_bin, :n_samples]
    required_params = [:κ, :log_A1, :E_c1, :γ, :log_A2, :E_c2, :Emin, :J1, :J2, :JE1, :JE2]
    required_suffixes = ["_median", "_q25", "_q75"]

    missing_cols = Symbol[]
    for col in required_base_cols
        if !hasproperty(df, col)
            push!(missing_cols, col)
        end
    end

    for param in required_params
        for suffix in required_suffixes
            col = Symbol(string(param) * suffix)
            if !hasproperty(df, col)
                push!(missing_cols, col)
            end
        end
    end

    if !isempty(missing_cols)
        error("DataFrame is missing required columns: $(join(missing_cols, ", "))")
    end

    # Create spatial bins named tuple
    spatial_bins = (
        mlat=collect(Float64, mlat_bins),
        mlt=collect(Float64, mlt_bins),
        ae=collect(String, ae_bins)
    )

    return EmpiricalFluxModel(
        df,
        spatial_bins,
        Tuple(Float64.(energy_range)),
        description,
        version
    )
end
