module DmspElfinConjunction

using Dates
using DataFrames
using DimensionalData

export get_elfin_flux_by_mlat, integrate_diff_flux
export get_flux_by_mlat  # re-export from this module

using SPEDAS
using SpacePhysicsMakie: set_if_valid!
using GeoCotrans

include("AACGM.jl")
include("fit.jl")
include("plot.jl")
include("makie.jl")

ntime(x) = size(x, 1)

function integrate_diff_flux(flux)
    Es = flux.dims[2].val
    return map(eachslice(flux, dims = 1)) do slc
        slc[1] * (Es[2] - Es[1]) +
            sum(slc[2:(end - 1)] .* (Es[3:end] .- Es[1:(end - 2)]) ./ 2) +
            slc[end] * (Es[end] - Es[end - 1])
    end
end

function get_flux_by_mlat(flux, mlat)
    # TODO: Improve the MLAT resolution by interpolating to 1 second first
    # Define MLAT bins (0.5° resolution)
    mlat_min = floor(minimum(mlat) * 2) / 2  # Round down to nearest 0.5
    mlat_max = ceil(maximum(mlat) * 2) / 2   # Round up to nearest 0.5
    mlat_bins = mlat_min:0.5:(mlat_max - 0.5)

    times = mlat.dims[1]

    res = map(mlat_bins) do bin
        idxs = findall(x -> bin <= x < bin + 0.5, mlat)
        if isempty(idxs)
            # Return missing or default values if no indices found
            return (missing, missing, missing, 0, 0)
        end
        min_idx = first(idxs)
        max_idx = min(last(idxs) + 1, length(times))
        mlat_t0 = times[min_idx]
        mlat_t1 = times[max_idx]
        flux_by_mlat = tview(flux, mlat_t0, mlat_t1)
        mean_flux = tmean(flux_by_mlat)
        mlat_t0, mlat_t1, mean_flux, ntime(flux_by_mlat), count(!isnan, mean_flux)
    end

    return DataFrame(;
        mlat = mlat_bins,
        mlat_t0 = getindex.(res, 1),
        mlat_t1 = getindex.(res, 2),
        flux = getindex.(res, 3),
        n_time = getindex.(res, 4),
        nnan_count = getindex.(res, 5)
    )
end


function get_flux_by_mlat(flux, mlat, timerange)
    flux_subset = tview(flux, timerange)
    mlat_subset = tview(mlat, timerange)
    return get_flux_by_mlat(flux_subset, mlat_subset)
end

# higher resolution of MLAT
function get_elfin_flux_by_mlat(flux, pos_gei, timerange)
    pos_aacgm = gei2aacgm(tview(pos_gei, timerange...))
    mlat = pos_aacgm.mlat
    return get_flux_by_mlat(flux, mlat, timerange)
end


function get_dmsp_flux_by_mlat(timerange, id)
    dmsp_flux = DMSP.load(timerange, id, "el_d_flux")
    dmsp_aacgm = DMSP.aacgm(timerange, id)
    dmsp_mlat_highres = getindex.(dmsp_aacgm, 1)
    return get_flux_by_mlat(dmsp_flux, dmsp_mlat_highres, timerange)
end

end
