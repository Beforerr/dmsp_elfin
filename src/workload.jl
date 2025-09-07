using Statistics: mean
using Speasy
import TimeseriesUtilities
using TimeseriesUtilities: tresample, tview
using PartialFunctions

"""
    gei2mlt_mlat(gei)
    gei2mlt_mlat(gei, timerange)

Convert ELFIN GEI coordinates to MLT and MLAT.

# Arguments
- `gei`: GEI coordinate data from ELFIN
- `timerange`: Optional time range to subset the data

# Returns
- Named tuple with `mlt` and `mlat` data with appropriate metadata
"""
gei2mlt_mlat(gei) = begin
    geo = gei2geo(gei)
    mlt = get_mlt(geo)
    aacgm = geo2aacgm(geo)
    mlat = aacgm.mlat
    return (
        mlt = setmeta(mlt, :ylabel => "MLT", :label => "ELFIN"),
        mlat = setmeta(mlat, :label => "ELFIN", :ylabel => "MLAT"),
    )
end

gei2mlt_mlat(gei, timerange) = gei2mlt_mlat(tview(gei, timerange...))

"""
    get_mlt_mlat(timerange, id)

Get MLT and MLAT data for DMSP `id` and `timerange` from SSC.

# Notes
- Data is fetched from SSC with 1-minute resolution and upsampled to 1-second
"""
get_mlt_mlat(timerange, id) = begin
    # the time resolution from SSC is 1 minute, floored/ceiled to minute boundaries
    t0 = floor(timerange[1], Minute)
    t1 = ceil(timerange[2], Minute)
    _geo = Speasy.get_data("ssc/dmspf$id/geo", t0, t1 + Second(1)) |> DimArray
    # Upsample to 1 second
    geo = tinterp(_geo, time_grid((t0, t1), Second(1)))
    mlt = get_mlt(geo)
    aacgm = geo2aacgm(geo)
    mlat = aacgm.mlat
    return mlt, mlat
end

# Remove the first channel of the flux if it is NaN
remove_first_channel(flux) = isnan(flux[1]) ? flux[2:end] : flux

"""
    workload(trange, ids; Δt=Minute(10), elx_flux=nothing, elx_gei=nothing, elx_probe="b", Δmlt_max=1, flux_threshold=200)

Process ELFIN and DMSP data for conjunction analysis and spectral fitting.

# Arguments
- `trange`: Time range tuple (start, end)
- `ids`: Vector of DMSP satellite IDs to process
- `Δt`: Time extension for DMSP data around the time range (default: 10 minutes)
- `elx_flux`: Pre-loaded ELFIN flux data (optional)
- `elx_gei`: Pre-loaded ELFIN GEI data (optional)
- `elx_probe`: ELFIN probe ("a" or "b", default: "b")
- `Δmlt_max`: Maximum MLT difference for conjunction (default: 1 hour)
- `flux_threshold`: Minimum flux threshold for quality filtering (default: 200)

# Returns
- `df`: Combined DataFrame with fitted spectral parameters for all valid conjunctions
- `elfin`: Named tuple with ELFIN data (flux, mlt, mlat)
- `dmsps`: Vector of named tuples with DMSP data for each satellite

# Processing Steps
1. Load/use ELFIN precipitating flux and position data
2. For each DMSP satellite:
   - Load flux data with time extension Δt
   - Calculate MLT/MLAT from orbital data
   - Bin data by MLAT (0.5° resolution)
   - Join ELFIN and DMSP data by MLAT
   - Filter by MLT difference and flux quality (≥3 channels > flux_threshold)
   - Fit two-step spectral models
3. Return combined results for analysis
"""
function workload(trange, ids; Δt = Minute(10), elx_flux = nothing, elx_gei = nothing, elx_probe = "b", Δmlt_max = 1, flux_threshold = 200, verbose = false, kw...)

    elx_mlt_trange = extend(trange, Second(1)) # Extend a bit to cover flux time range

    elx_flux = @something elx_flux ELFIN.precipitating_flux(trange, elx_probe)
    elx_gei = @something elx_gei ELFIN.gei(elx_mlt_trange, elx_probe)

    elx_pflux = tview(elx_flux.prec, trange)
    elx_mlt, elx_mlat = gei2mlt_mlat(tview(elx_gei, elx_mlt_trange))
    elx_df = get_flux_by_mlat(elx_pflux, elx_mlat)

    elfin = (; flux = elx_pflux, mlt = elx_mlt, mlat = elx_mlat)

    results = map(ids) do id
        dmsp_range = extend(trange, Δt)
        f = DMSP.flux(dmsp_range, id)
        mlt, mlat = get_mlt_mlat(dmsp_range, id)
        if isnothing(f)
            @warn "Skipping DMSP ID $id due to missing flux data for range $dmsp_range"
            return (; df = DataFrame(), dmsp = nothing)
        end
        dmsp = (; id, flux = f, mlt, mlat)
        dmsp_df = get_flux_by_mlat(f, mlat)

        df = @chain begin
            leftjoin(elx_df, dmsp_df, on = :mlat, renamecols = "_elx" => "_dmsp")
            sort!(:mlat)
            dropmissing!()
            @rtransform! @astable begin
                :success = false
                :dmsp_mlt = mean(tview(mlt, :mlat_t0_dmsp, :mlat_t1_dmsp))
                :elfin_mlt = mean(tview(elx_mlt, :mlat_t0_elx, :mlat_t1_elx))
                :Δmlt = mlt_dist(:dmsp_mlt, :elfin_mlt)
                :id = id
            end
            # Notes: there is a subset of orbits with no working first channel. So we remove the first channel.
            @rtransform! :flux_elx = remove_first_channel(:flux_elx)
            # require at least having 3 lowest channels on elfin should have flux larger than 200
            @aside @info "There are $(nrow(_)) elements before filtering"
            @rsubset! begin
                :Δmlt < Δmlt_max
                # sum(:flux_elx[1:3] .> 200) >= 3
            end
            @aside @info "There are $(nrow(_)) elements after filtering"
            @rtransform! $AsTable = begin
                elx_flux = :flux_elx
                # remove all ELFIN measurements after the first (lowest energy) NaN
                nan_idx = findfirst(isnan, elx_flux)
                elx_flux = isnothing(nan_idx) ? elx_flux : elx_flux[1:(nan_idx - 1)]
                fit_two_flux(:flux_dmsp, elx_flux; verbose, kw...)
            end
            @rsubset! :success
            @aside @info "There are $(nrow(_)) successful fits"
        end
        (; df, dmsp)
    end

    # Filter out failed results and combine successful ones
    successful_results = filter(x -> !isempty(x.df), results)

    if isempty(successful_results)
        @warn "No successful DMSP data processing for any ID in $ids"
        return DataFrame(), elfin, []
    end

    df = mapreduce(first, vcat $ (; cols = :union), successful_results)
    dmsps = map(last, successful_results)
    return df, elfin, dmsps
end


"""
    model_flux(model, mlat)

Generate model flux data on a continuous MLAT grid from min to max MLAT with gaps filled by NaN. 
Returns `DimArray`: 2D array with dimensions (MLAT, Energy) containing model flux values

# Example
If mlat = [-74.0, -72.0, -71.0], the function creates a grid:
[-74.0, -73.5, -73.0, -72.5, -72.0, -71.5, -71.0]
where -73.5, -73.0, -72.5, -71.5 are filled with NaN.
"""
model_flux(models, mlat; δmlat = 0.5) = begin
    Es = default_energies()

    mlat_min, mlat_max = extrema(mlat)
    full_mlat = mlat_min:δmlat:mlat_max  # 0.5° bins as edges

    # Create data matrix with NaNs for missing MLATs
    data = fill(NaN, length(full_mlat), length(Es))

    # Fill in data for available MLATs
    for (i, m) in enumerate(full_mlat)
        # Find closest MLAT bin edge
        idx = findfirst(==(m), mlat)
        if !isnothing(idx)
            model = models[idx]
            data[i, :] .= model.(Es)
        end
    end

    dims = (X(full_mlat), Y(Es))
    return DimArray(data, dims)
end

axis(f) = f.axis

"""
    plot_flux_analysis(df, elfin, dmsps; figsize=(1200, 1000))

Plot flux analysis with ELFIN and multiple DMSP data, including modeled fluxes and parameter variations.

- `elfin`: Named tuple with ELFIN data (flux, mlt, mlat)  
- `dmsps`: Vector of named tuples with DMSP data (flux, mlt, mlat)

# Returns
- Figure with flux plots, modeled fluxes
"""
function plot_flux_analysis(f, df, elfin, dmsps; mlats = nothing)
    colorrange = COLORRANGE[]

    # Plot MLAT - MLT
    mlt_ax = (; xlabel = "MLAT", ylabel = "MLT")
    mlt_axplot = lines(f[1, 1:2], elfin.mlat.data, elfin.mlt.data; label = "ELFIN", axis = mlt_ax)

    # Plot ELFIN flux
    replace!(elfin.flux, 0 => NaN)
    p_elfin = plot_flux_by_mlat(f[2, 1:2], elfin.flux, elfin.mlat; colorrange)
    xlims = extrema(elfin.mlat)

    # Plot DMSP fluxes
    idx = 3
    dmsp_plots = mapreduce(vcat, enumerate(dmsps)) do (i, dmsp)
        @unpack id, mlat, flux = dmsp
        gdf = @subset(df, :id .== id)
        flux = replace!(flux, 0 => NaN)
        p1 = plot_flux_by_mlat(f[idx, 1:2], flux, mlat; colorrange)
        xlims!(p1.axis, xlims)
        modeled_fluxes = model_flux(gdf.model, gdf.mlat)

        p2 = plot_flux_by_mlat(f[idx + 1, 1:2], modeled_fluxes; colorrange)
        idx += 2

        # Plot MLAT - MLT
        lines!(mlt_axplot.axis, mlat.data, dmsp.mlt.data; label = "DMSP $id")

        if !isnothing(mlats)
            sdf = @rsubset(gdf, :mlat ∈ mlats)
            if nrow(sdf) > 0
                plot_spectra(f[idx, 1:2], sdf)
                idx += 1
            end
        end

        hidexdecorations!(p1.axis; grid = false)

        [p1, p2]
    end
    flux_plots = [p_elfin, dmsp_plots...]

    axs = axis.(flux_plots)
    linkxaxes!(axs...)

    Colorbar(f[2:end, 3], p_elfin.plot; label = 𝒀.nflux)
    if !isnothing(mlats)
        vlines!.(axs, Ref(mlats))
    end


    all_plots = [flux_plots..., mlt_axplot]
    # Set x-axis hide decorations
    map(all_plots) do p
        xlims!(p.axis, xlims)
    end

    map((mlt_axplot, p_elfin)) do p
        hidexdecorations!(p.axis; grid = false)
    end

    axislegend(mlt_axplot.axis, position = :lb)

    # # Plot parameter variations if fits are successful
    # if !isempty(successful_fits)
    #     param_start_col = 4
    #     @with successful_fits begin
    #         plot_parameters_variation(f[1:4+length(dmsps)-1, param_start_col:param_start_col+1],
    #                                 :mlat, :model, :n_points; scores = :score)
    #     end
    # end

    return all_plots
end

function plot_flux_analysis(args...; kw...)
    f = Figure(; size = (1200, 1000))
    plots = plot_flux_analysis(f, args...; kw...)
    return f, plots
end


function plot_elfin_dmsp(timerange, ids; elx_flux, elx_gei, Δt = Minute(10))
    ## Part 1: Fluxes by Time
    dmsp_fluxs = map(ids) do id
        DMSP.flux(extend(timerange, Δt), id)
    end

    elx_flux = tview(elx_flux, timerange)

    foreach((elx_flux, dmsp_fluxs...)) do flux
        set_if_valid!(
            flux.metadata,
            :yscale => log10, :colorscale => log10, :colorrange => COLORRANGE[]
        )
        replace!(flux, 0 => NaN)
    end

    elx_mlt, elx_mlat = gei2mlt_mlat(tview(elx_gei, extend(timerange, Second(1)))) # extend so that interpolation works

    xlims = extrema(elx_mlat)

    dmsp_mlt_mlats = map(ids) do id
        mlt, mlat = get_mlt_mlat(extend(timerange, Δt + Second(2)), id)
        setmeta(mlt, :label => "DMSP $id", :ylabel => "MLT"), setmeta(mlat, :label => "DMSP $id", :ylabel => "MLAT")
    end
    mlts = [elx_mlt, first.(dmsp_mlt_mlats)...]
    mlats = [elx_mlat, last.(dmsp_mlt_mlats)...]

    tvars = (elx_flux, dmsp_fluxs..., mlts, mlats)
    faxs = SpacePhysicsMakie.tplot(tvars)

    ## Part 2: Fluxes by MLAT
    f = faxs.figure
    elx_axplot = plot_flux_by_mlat(f[1, 2], elx_flux, elx_mlat)
    hidexdecorations!(elx_axplot.axis)
    for i in eachindex(dmsp_fluxs)
        axplot = plot_flux_by_mlat(f[1 + i, 2], dmsp_fluxs[i], mlats[1 + i])
        ax = axplot.axis
        xlims!(ax, xlims)
        hidexdecorations!(ax)
    end

    markersize = 2
    ax = Axis(f[length(dmsp_fluxs) + 2, 2])
    scatterlines!(ax, elx_mlat.data, elx_mlt.data; label = "ELFIN", markersize)
    foreach(dmsp_mlt_mlats) do (mlt, mlat)
        scatterlines!(ax, mlat.data, mlt.data; label = "$(mlat.metadata[:label])", markersize)
    end
    axislegend(ax)
    xlims!(ax, xlims)

    return faxs
end
