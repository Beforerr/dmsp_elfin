using Statistics: mean
using Speasy
import TimeseriesUtilities
using TimeseriesUtilities: tresample, tview, tsort, find_continuous_timeranges
using PartialFunctions
using DrWatson
using Logging
using LoggingExtras: TeeLogger
using DmspElfinConjunction: get_trange_by_mlat
import DmspElfinConjunction as DE
import ELFINData as ELFIN
using DimensionalData: NoMetadata
import DimensionalData as DD
using SpectralModels: fit_two_flux

include("./conjugation.jl")
include("./mechanisms.jl")

const FLUX_THRESHOLD = 150

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

get_geo(id, t0, t1) = DimArray(Speasy.get_data("ssc/dmspf$id/geo", t0, t1))

"""
    get_mlt_mlat(id, timerange)

Get MLT and MLAT data for DMSP `id` and `timerange` from SSC.

# Notes
- Data is fetched from SSC with 1-minute resolution and upsampled to 1-second
"""
get_mlt_mlat(id, t0, t1) = begin
    # the time resolution from SSC is 1 minute, floored/ceiled to minute boundaries
    t0 = floor(DateTime(t0), Minute)
    t1 = ceil(DateTime(t1), Minute)
    dt = Second(1)
    _geo = get_geo(id, t0, t1 + dt)
    # Upsample to 1 second
    geo = tinterp(_geo, t0:dt:t1)
    mlt = get_mlt(geo)
    aacgm = geo2aacgm(geo)
    mlat = aacgm.mlat
    return mlt, mlat
end

get_mlt_mlat(id, trange) = get_mlt_mlat(id, trange[1], trange[2])

"""
    sanitize_dmsp_flux(flux)

1. Remove the last channel
"""
function sanitize_dmsp_flux(flux)
    return @views flux[1:(end - 1)]
end
# The difference (uncertainty) in kappa between just removing the last channel and doing nothing is about 0.2-2 (usually 0.5).


function process_dmsp_conjunction(id, trange, elx_df; Δmlt_max = 1, kw...)
    mlt, mlat = get_mlt_mlat(id, trange)
    dmsp_df = DE.get_trange_by_mlat(mlat)
    df = @chain begin
        leftjoin(elx_df, dmsp_df, on = :mlat, renamecols = "_elx" => "_dmsp")
        dropmissing!()
        @rtransform! @astable begin
            :success = false
            :mlt_dmsp = local_mlt_mean(tview(mlt, :trange_dmsp))
            :Δmlt = mlt_dist(:mlt_dmsp, :mlt_elx)
            :id = id
        end
        # require at least having 3 lowest channels on elfin should have flux larger than 200
        @aside @info "There are $(nrow(_)) elements before filtering"
        @rsubset! :Δmlt < Δmlt_max
    end

    @info "There are $(nrow(df)) elements after filtering"
    nrow(df) == 0 && return DataFrame()

    flux = DMSP.flux(trange, id)
    if isnothing(flux) || isempty(flux)
        @warn "Skipping DMSP ID $id due to missing or empty flux data for range $trange"
        return DataFrame()
    end
    flux = tsort(flux)
    df = @chain df begin
        @rtransform! @astable begin
            flux_by_mlat = tview(flux, :trange_dmsp)
            :flux_dmsp = tmean(flux_by_mlat)
            :n_time_dmsp = DE.ntime(flux_by_mlat)
        end
        @rtransform! $AsTable = fit_two_flux(sanitize_dmsp_flux(:flux_dmsp), :flux_elx; threshold = FLUX_THRESHOLD, kw...)
        @rsubset! :success
    end
    @info "There are $(nrow(df)) successful fits"
    return df
end

"""
1. Remove the first channel of the flux if it is NaN or lower than second channel (this is a problem with the first channel only for some events)
2. Remove all ELFIN measurements after the first (lowest energy) NaN
"""
function sanitize_elfin_flux(flux)
    _flux = (isnan(flux[1]) || flux[1] < flux[2]) ? flux[2:end] : flux
    nan_idx = findfirst(isnan, _flux)
    return isnothing(nan_idx) ? _flux : _flux[1:(nan_idx - 1)]
end

"""
    workload(trange, ids; Δt=Minute(10), elx_flux=nothing, elx_gei=nothing, elx_probe="b", Δmlt_max=1, flux_threshold=200)

Process ELFIN and DMSP data for conjunction analysis and spectral fitting.

# Arguments
- `trange`: Time range tuple (start, end)
- `ids`: Vector of DMSP satellite IDs to process
- `Δt`: Time extension for DMSP data around the time range (default: 10 minutes)
- `Δmlt_max`: Maximum MLT difference for conjunction (default: 1 hour)
- `flux_threshold`: Minimum flux threshold for quality filtering (default: 200)

`elx_flux` and `elx_gei` are optional arguments for pre-loaded ELFIN data. If not provided, they will be loaded automatically for ELFIN probe `elx_probe` ("a" or "b").

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
   - Filter by MLT difference and flux quality (≥3 channels > FLUX_THRESHOLD)
   - Fit two-step spectral models
3. Return combined results for analysis

See also: `fit_two_flux`
"""
function workload(trange, ids, elx_flux, elx_gei; Δt = Minute(10), kw...)
    elx_mlt_trange = extend(trange, Second(1)) # Extend a bit to cover flux time range
    elx_mlt, elx_mlat = gei2mlt_mlat(tview(elx_gei, elx_mlt_trange))

    elx_df = @chain get_trange_by_mlat(elx_mlat) begin
        @rtransform! begin
            :mlt = local_mlt_mean(tview(elx_mlt, :trange))
            :flux = tmean(tview(elx_flux.prec, :trange)) |> sanitize_elfin_flux
        end
        @rsubset! length(:flux) >= 3 && all(view(:flux, 1:3) .> FLUX_THRESHOLD)
    end

    results = map(ids) do id
        process_dmsp_conjunction(id, extend(trange, Δt), elx_df; kw...)
    end
    successful_results = filter(!isempty, results)

    if isempty(successful_results)
        @warn "No successful DMSP data processing for any ID in $ids"
        return DataFrame()
    end
    df = reduce(vcat, successful_results)
    classify_precipitation!(df, elx_flux, elx_mlt; trange = :trange_elx, kw...)
    return df
end

workload((trange, ids); kw...) = workload(trange, ids; kw...)

tspan(timerange) = timerange[2] - timerange[1]

function delete_metadata!(df)
    return @rtransform! df begin
        :flux_elx = rebuild(:flux_elx, metadata = NoMetadata())
        :flux_dmsp = rebuild(:flux_dmsp, metadata = NoMetadata())
    end
end

function add_metadata!(df, conf)
    return isempty(df) ? df : begin
            @rtransform! df begin
                :flux_elx = rebuild(:flux_elx, metadata = conf["elfin_metadata"])
                :flux_dmsp = rebuild(:flux_dmsp, metadata = conf["dmsps_metadata"])
            end
        end
end

"""
    produce(trange, probe, ids; Δt = Minute(10), Δmlt_max = 1, kw...)

Process ELFIN and DMSP statistic results for analysis.

See also: [`workload`](@ref)
"""
function produce(trange, probe, ids; Δt = Minute(10), Δmlt_max = 1, force = false, kw...)
    conf = @strdict trange probe ids
    filename = "ELFIN=$(probe)_t0=$(trange[1])_t1=$(trange[2])_ids=$(join(ids, "-"))"
    return produce_or_load(conf, datadir(); filename, force) do c
        trange, probe = c["trange"], c["probe"]
        elx_gei = ELFIN.gei(trange, probe)
        elx_flux = tsort(permutedims(ELFIN.epd(trange, probe)))
        continuous_ranges = find_continuous_timeranges(elx_flux, Second(60))
        # filter very short ranges
        filter!(x -> tspan(x) > Second(5), continuous_ranges)

        valid_ranges_with_ids = find_matched_mlat_conditions(
            continuous_ranges, elx_gei; ids, Δt, Δmlt_max
        )
        dfs = map(valid_ranges_with_ids) do (trange, ids)
            df = workload(trange, ids, elx_flux, elx_gei; Δt)
            df.maxAE .= maxAE(trange)[1]
            return df
        end
        filter!(!isempty, dfs)
        df = isempty(dfs) ? DataFrame() : reduce(vcat, dfs)
        conf["df"] = df
        isempty(df) ? conf : begin
                delete_metadata!(df)
                conf["elfin_metadata"] = DD.metadata(elx_flux.prec)
                conf["dmsps_metadata"] = DD.metadata(df.flux_dmsp[1])
                conf
            end
    end
end

function produce(tranges::AbstractVector, probe, ids; file = "logfile.txt", kw...)
    return open(file, "a+") do io
        file_logger = SimpleLogger(io)
        tee_logger = TeeLogger(current_logger(), file_logger)

        with_logger(tee_logger) do
            dfs = map(tranges) do trange
                conf = produce(trange, probe, ids; kw...)[1]
                add_metadata!(conf["df"], conf)
            end
            reduce(vcat, filter(!isempty, dfs))
        end
    end
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
model_flux(models, mlats; δmlat = 0.5) = begin
    Es = default_energies()

    mlat_min, mlat_max = extrema(mlats)
    full_mlat = mlat_min:δmlat:mlat_max  # 0.5° bins as edges

    # Create data matrix with NaNs for missing MLATs
    data = fill(NaN, length(full_mlat), length(Es))

    # Fill in data for available MLATs
    for (i, m) in enumerate(full_mlat)
        # Find closest MLAT bin edge
        idx = findfirst(==(m), mlats)
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
    elx_flux = replace(elfin.flux, 0 => NaN)
    p_elfin = plot_flux_by_mlat(f[2, 1:2], elx_flux, elfin.mlat; colorrange)
    xlims = extrema(elfin.mlat)

    # Plot DMSP fluxes
    idx = 3
    dmsp_plots = mapreduce(vcat, dmsps) do dmsp
        @unpack id, mlat = dmsp
        gdf = @subset(df, :id .== id)
        flux = replace(dmsp.flux, 0 => NaN)
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
    all_plots = [mlt_axplot, p_elfin, dmsp_plots...]

    axs = axis.(all_plots)
    linkxaxes!(axs...)

    Colorbar(f[2:end, 3], p_elfin.plot; label = 𝒀.nflux)
    !isnothing(mlats) && vlines!.(axs, Ref(mlats))
    # Set x-axis hide decorations
    xlims!.(axs, Ref(xlims))

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
        # replace!(flux, 0 => NaN)
    end

    elx_mlt, elx_mlat = gei2mlt_mlat(tview(elx_gei, extend(timerange, Second(1)))) # extend so that interpolation works

    xlims = extrema(elx_mlat)

    dmsp_mlt_mlats = map(ids) do id
        mlt, mlat = get_mlt_mlat(id, extend(timerange, Δt + Second(2)))
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
