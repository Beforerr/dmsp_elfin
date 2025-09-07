module DmspElfinConjunction

using Dates
using DataFrames
using DimensionalData
import Speasy

export get_elfin_flux_by_mlat, integrate_diff_flux
export get_flux_by_mlat  # re-export from this module

using SPEDAS
using GeoCotrans
using Printf
export energies
export dist, mlt_dist
export maxAE
export extend

include("utils.jl")
include("AACGM.jl")
include("fit.jl")
include("makie.jl")

ntime(x) = size(x, 1)
extend(timerange, Δt) = (timerange[1] - Δt, timerange[2] + Δt)

# use maxAE for three hours before the event
function maxAE(trange, dt = Hour(3))
    pre_range = DateTime.(trange) .- dt
    ae = Speasy.get_data("cda/OMNI_HRO2_1MIN/AE_INDEX", pre_range...)
    return maximum(ae), ae
end

function integrate_diff_flux(flux)
    Es = flux.dims[2].val
    return map(eachslice(flux, dims = 1)) do slc
        slc[1] * (Es[2] - Es[1]) +
            sum(slc[2:(end - 1)] .* (Es[3:end] .- Es[1:(end - 2)]) ./ 2) +
            slc[end] * (Es[end] - Es[end - 1])
    end
end

function get_flux_by_mlat(flux, mlats; δmlat = 0.5, δt = Millisecond(1000))
    # TODO: Improve the MLAT resolution by interpolating to 1 second first
    # left closed, right open
    mlat_min = round(minimum(mlats) * 2, RoundNearestTiesUp) / 2  # Round down to nearest 0.5
    mlat_max = round(maximum(mlats) * 2, RoundNearestTiesUp) / 2   # Round up to nearest 0.5
    mlat_bins = mlat_min:δmlat:(mlat_max - δmlat)

    times = mlats.dims[1].val
    resolution(times) == δt || error("Time resolution does not match δt")

    res = map(mlat_bins) do bin
        idxs = findall(x -> bin - δmlat / 2 <= x < bin + δmlat / 2, mlats)
        if isempty(idxs)
            # Return missing or default values if no indices found
            return (missing, missing, missing, 0, 0)
        end
        min_idx = first(idxs)
        max_idx = min(last(idxs), length(times))
        mlat_t0 = times[min_idx] - δt / 2
        mlat_t1 = times[max_idx] + δt / 2
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
