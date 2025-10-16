"""
Model creation and parameter extraction functions.
"""

"""
    get_ae_bin(ae_value, ae_bins::Vector{String}) -> String

Find the appropriate AE bin for a given AE index value.

# Example
```julia
bins = ["[0, 100)", "[100, 300)", "[300, Inf)"]
get_ae_bin(150, bins)  # Returns "[100, 300)"
```
"""
function get_ae_bin(ae_value, ae_bins::Vector{String})
    for bin in ae_bins
        # Parse bin string like "[0, 100)" or "[300, Inf)"
        m = match(r"\[(\d+(?:\.\d+)?),\s*(\d+(?:\.\d+)?|Inf)\)", bin)
        if m !== nothing
            lower = parse(Float64, m.captures[1])
            upper = m.captures[2] == "Inf" ? Inf : parse(Float64, m.captures[2])
            if lower <= ae_value < upper
                return bin
            end
        end
    end
    # Return closest bin if exact match not found
    return ae_bins[end]  # Default to highest AE bin
end

"""
    interpolate_parameters(model::EmpiricalFluxModel, mlat::Real, mlt::Real, ae::Real;
                          stat::Symbol=:median) -> FluxParameters

Interpolate model parameters for arbitrary geophysical conditions.

Uses nearest-neighbor interpolation in (MLat, MLT, AE) space.

# Arguments
- `model`: The empirical flux model
- `mlat`: Magnetic latitude in degrees
- `mlt`: Magnetic local time in hours
- `ae`: AE index in nT
- `stat`: Which statistic to use (`:median`, `:q25`, `:q75`)

# Returns
`FluxParameters` object containing all spectral model parameters.
"""
function interpolate_parameters(
        model::EmpiricalFluxModel, mlat::Real, mlt::Real, ae::Real;
        stat::Symbol = :median
    )
    # Get appropriate AE bin
    ae_bin = get_ae_bin(ae, model.spatial_bins.ae)

    # Filter for AE bin
    df_ae = @rsubset(model.parameters, :ae_bin == ae_bin)

    if nrow(df_ae) == 0
        error("No data available for AE bin: $ae_bin")
    end

    # Find nearest neighbor in (MLat, MLT) space
    # Normalize coordinates for distance calculation
    mlat_norm = abs(mlat)
    mlt_norm = mod(mlt, 24.0)

    # Handle MLT wrapping (0h ≈ 24h)
    distances = map(eachrow(df_ae)) do row
        Δmlat = row.mlat_bin - mlat_norm
        Δmlt = min(
            abs(row.mlt_bin - mlt_norm),
            abs(row.mlt_bin - mlt_norm + 24),
            abs(row.mlt_bin - mlt_norm - 24)
        )
        sqrt(Δmlat^2 + (Δmlt * 2)^2)  # Weight MLT less than MLat
    end

    nearest_idx = argmin(distances)
    row = df_ae[nearest_idx, :]

    # Extract parameters with requested statistic
    suffix = stat == :median ? "_median" : (stat == :q25 ? "_q25" : "_q75")

    return FluxParameters(
        row[Symbol("log_A1" * suffix)],
        row[Symbol("E_c1" * suffix)],
        row[Symbol("γ" * suffix)],
        row[Symbol("log_A2" * suffix)],
        row[Symbol("E_c2" * suffix)],
        row[Symbol("κ" * suffix)],
        row[Symbol("Emin" * suffix)],
        row[Symbol("J1" * suffix)],
        row[Symbol("J2" * suffix)],
        row[Symbol("JE1" * suffix)],
        row[Symbol("JE2" * suffix)],
        row.mlat_bin,
        row.mlt_bin,
        ae,
        row.n_samples
    )
end

"""
    flux_parameters(model::EmpiricalFluxModel; mlat, mlt, ae, stat=:median) -> FluxParameters

Get flux model parameters for specified conditions.

Convenience wrapper around `interpolate_parameters`.
"""
flux_parameters(model::EmpiricalFluxModel; mlat::Real, mlt::Real, ae::Real, stat::Symbol = :median) =
    interpolate_parameters(model, mlat, mlt, ae; stat = stat)
