"""
Script to create PrecipitatingFluxModels from statistical analysis results.

This script processes the statistical analysis results from the DMSP-ELFIN
conjunction study and creates a serialized empirical flux model for end users.
"""

using DataFrames, DataFramesMeta
using DmspElfinConjunction
using SpectralModels: A, κ, E_c, n_flux, e_flux
using Statistics
using CategoricalArrays

# Load the PrecipitatingFluxModels package
using PrecipitatingFluxModels

PrecipitatingFluxModels.save("default_model.jld2", model)
PrecipitatingFluxModels.load()

"""
    prepare_model_parameters(sdf::DataFrame) -> DataFrame

Prepare statistical parameters DataFrame for the flux model.

Transforms the analysis results into the format required by PrecipitatingFluxModels.
"""
function prepare_model_parameters(gdf::DataFrame)
    # Compute statistics for all parameters
    params_df = combine(
        gdf,
        nrow => :n_samples,
        # Kappa parameter
        :κ => median => :κ_median,
        :κ => (x -> quantile(x, 0.25)) => :κ_q25,
        :κ => (x -> quantile(x, 0.75)) => :κ_q75,
        # ExpPow model (model1)
        :log_A1 => median => :log_A1_median,
        :log_A1 => (x -> quantile(x, 0.25)) => :log_A1_q25,
        :log_A1 => (x -> quantile(x, 0.75)) => :log_A1_q75,
        :E_c1 => median => :E_c1_median,
        :E_c1 => (x -> quantile(x, 0.25)) => :E_c1_q25,
        :E_c1 => (x -> quantile(x, 0.75)) => :E_c1_q75,
        :γ => median => :γ_median,
        :γ => (x -> quantile(x, 0.25)) => :γ_q25,
        :γ => (x -> quantile(x, 0.75)) => :γ_q75,
        # Kappa model (model2)
        :log_A2 => median => :log_A2_median,
        :log_A2 => (x -> quantile(x, 0.25)) => :log_A2_q25,
        :log_A2 => (x -> quantile(x, 0.75)) => :log_A2_q75,
        :E_c2 => median => :E_c2_median,
        :E_c2 => (x -> quantile(x, 0.25)) => :E_c2_q25,
        :E_c2 => (x -> quantile(x, 0.75)) => :E_c2_q75,
        # Transition energy
        :Emin => median => :Emin_median,
        :Emin => (x -> quantile(x, 0.25)) => :Emin_q25,
        :Emin => (x -> quantile(x, 0.75)) => :Emin_q75,
        # Number fluxes
        :J1 => median => :J1_median,
        :J1 => (x -> quantile(x, 0.25)) => :J1_q25,
        :J1 => (x -> quantile(x, 0.75)) => :J1_q75,
        :J2 => median => :J2_median,
        :J2 => (x -> quantile(x, 0.25)) => :J2_q25,
        :J2 => (x -> quantile(x, 0.75)) => :J2_q75,
        # Energy fluxes
        :JE1 => median => :JE1_median,
        :JE1 => (x -> quantile(x, 0.25)) => :JE1_q25,
        :JE1 => (x -> quantile(x, 0.75)) => :JE1_q75,
        :JE2 => median => :JE2_median,
        :JE2 => (x -> quantile(x, 0.25)) => :JE2_q25,
        :JE2 => (x -> quantile(x, 0.75)) => :JE2_q75,
    ) |> dropmissing!

    # Convert bin labels to numeric values
    @transform!(
        params_df,
        :mlat_bin = _mlat.(:mlat_bin),
        :mlt_bin = parse.(Int, :mlt_bin),
        :ae_bin = string.(:maxAE_bin)
    )

    # Remove the temporary maxAE_bin column
    select!(params_df, Not(:maxAE_bin))

    return params_df
end

"""
    main()

Main function to create and save the flux model.
"""
function main()
    @info "Loading statistical analysis data..."

    # Load the analysis results (adjust path as needed)
    # This assumes you've run the statistics analysis and have `sdf` available
    include(joinpath(@__DIR__, "..", "notebooks", "stats.qmd"))

    @info "Found $(nrow(sdf)) observations"

    @info "Preparing model parameters..."
    params_df = prepare_model_parameters(sdf)

    @info "Created parameter table with $(nrow(params_df)) spatial bins"

    # Define spatial bins
    mlat_bins = sort(unique(params_df.mlat_bin))
    mlt_bins = sort(unique(params_df.mlt_bin))
    ae_bins = sort(unique(params_df.ae_bin))

    @info "Spatial coverage:"
    @info "  MLat: $(length(mlat_bins)) bins ($(minimum(mlat_bins))° - $(maximum(mlat_bins))°)"
    @info "  MLT: $(length(mlt_bins)) bins ($(minimum(mlt_bins))h - $(maximum(mlt_bins))h)"
    @info "  AE: $(length(ae_bins)) bins: $(join(ae_bins, ", "))"

    # Create model
    @info "Creating empirical flux model..."
    model = create_model_from_stats(
        params_df;
        mlat_bins = mlat_bins,
        mlt_bins = mlt_bins,
        ae_bins = ae_bins,
        energy_range = (0.03, 1000.0),
        description = "DMSP-ELFIN precipitating electron flux model (2020-2022)",
        version = "1.0.0"
    )

    return model
end
