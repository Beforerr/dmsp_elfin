module DmspElfinConjunction

using Dates
using HDF5
using JLD2
using DataFrames, DataFramesMeta
using DimensionalData

export read2dimarray, read2dimstack
export get_flux_by_mlat, get_dmsp_flux_by_mlat, get_elfin_flux_by_mlat

using SPEDAS
using DataInterpolations: ExtrapolationType
using GeoCotrans

include("DMSP.jl")
include("hdf5.jl")
include("AACGM.jl")
include("fit.jl")

ntime(x) = size(x, 1)

function get_flux_by_mlat(flux, mlat, timerange)
    # Improve the MLAT resolution by interpolating to 1 second first
    flux_subset = tview(flux, timerange...)
    mlat_subset = tview(mlat, timerange...)

    # Define MLAT bins (0.5° resolution)
    mlat_min = floor(minimum(mlat_subset) * 2) / 2  # Round down to nearest 0.5
    mlat_max = ceil(maximum(mlat_subset) * 2) / 2   # Round up to nearest 0.5
    mlat_bins = mlat_min:0.5:mlat_max-0.5

    times = mlat_subset.dims[1]

    res = map(mlat_bins) do bin
        idxs = findall(x -> bin <= x < bin + 0.5, mlat_subset)
        min_idx = first(idxs)
        max_idx = min(last(idxs) + 1, length(times))
        mlat_t0 = times[min_idx]
        mlat_t1 = times[max_idx]
        flux_by_mlat = tview(flux_subset, mlat_t0, mlat_t1)
        mean_flux = tmean(flux_by_mlat)
        mlat_t0, mlat_t1, mean_flux, ntime(flux_by_mlat), count(!isnan, mean_flux)
    end

    return DataFrame(;
        mlat=mlat_bins,
        mlat_t0=getindex.(res, 1),
        mlat_t1=getindex.(res, 2),
        flux=getindex.(res, 3),
        n_time=getindex.(res, 4),
        nnan_count=getindex.(res, 5)
    )
end

# higher resolution of MLAT
function get_elfin_flux_by_mlat(flux, pos_gei, timerange)
    pos_aacgm = gei2aacgm(tview(pos_gei, timerange...))
    mlat = pos_aacgm.mlat
    get_flux_by_mlat(flux, mlat, timerange)
end


function get_dmsp_flux_by_mlat(timerange, id)
    dmsp_flux = dmsp_load(timerange, id, "el_d_flux")
    dmsp_aacgm = dmsp_get_aacgm(timerange, id)
    dmsp_mlat_highres = getindex.(dmsp_aacgm, 1)
    get_flux_by_mlat(dmsp_flux, dmsp_mlat_highres, timerange)
end


end
