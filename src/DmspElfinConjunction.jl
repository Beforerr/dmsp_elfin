module DmspElfinConjunction

using Dates
using DataFrames, DataFramesMeta
using DimensionalData
using Dictionaries
import Speasy
using NaNStatistics

export integrate_diff_flux
export get_flux_by_mlat

using SPEDAS
using GeoCotrans
using Printf
export energies
export dist, mlt_dist, local_mlt_mean
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
    return isempty(ae) ? missing : maximum(ae), ae
end

function integrate_diff_flux(flux)
    Es = flux.dims[2].val
    return map(eachslice(flux, dims = 1)) do slc
        slc[1] * (Es[2] - Es[1]) +
            sum(slc[2:(end - 1)] .* (Es[3:end] .- Es[1:(end - 2)]) ./ 2) +
            slc[end] * (Es[end] - Es[end - 1])
    end
end

"""
    bin_mlat_times(mlats; δmlat = 0.5, δt = Millisecond(1000))

Divide `mlats` into bins with bin size `δmlat` and return time ranges for each bin with time resolution `δt`.

Returns a dictionary mapping `mlat_bin => (t0, t1)` where t0 and t1 are the time bounds for that MLAT bin.
"""
function bin_mlat_times(mlats; δmlat = 0.5, δt = Millisecond(1000))
    # Create MLAT bins - left closed, right open
    mlat_min = round(nanminimum(mlats) * 2, RoundNearestTiesUp) / 2
    mlat_max = round(nanmaximum(mlats) * 2, RoundNearestTiesUp) / 2
    mlat_bins = mlat_min:δmlat:(mlat_max - δmlat)

    times = mlats.dims[1].val.data
    resolution(times) == δt || error("Time resolution does not match δt")

    # Create mapping of mlat_bin => (t0, t1)
    bin_times = map(mlat_bins) do bin
        idxs = findall(x -> bin - δmlat / 2 <= x < bin + δmlat / 2, mlats)
        N = length(idxs)
        if N == 0
            nothing
        elseif N == 1
            [(times[idxs[1]] - δt / 2, times[idxs[1]] + δt / 2)]
        else
            # Group consecutive indices into continuous segments
            segments = Tuple{Int, Int}[]
            current_start = idxs[1]

            for i in 2:N
                current_idx = idxs[i]
                previous_idx = idxs[i - 1]
                if current_idx != previous_idx + 1
                    push!(segments, (current_start, previous_idx))
                    current_start = idxs[i]
                end
            end
            # Don't forget the last segment
            push!(segments, (current_start, idxs[end]))

            # Convert index segments to time ranges
            map(segments) do (start_idx, end_idx)
                t0 = times[start_idx] - δt / 2
                t1 = times[end_idx] + δt / 2
                (t0, t1)
            end
        end
    end
    return filter(!isnothing, Dictionary(mlat_bins, bin_times))
end

function get_trange_by_mlat(mlats; kw...)
    bin_times = bin_mlat_times(mlats; kw...)
    # Expand each MLAT bin with multiple segments into separate rows
    return mapreduce(vcat, pairs(bin_times)) do (mlat_bin, time_ranges)
        DataFrame(;
            mlat = fill(mlat_bin, length(time_ranges)),
            trange = time_ranges,
        )
    end
end

function get_flux_by_mlat(flux, mlats; kw...)
    df = get_trange_by_mlat(mlats; kw...)
    return @rtransform! df @astable begin
        flux_by_mlat = tview(flux, :trange)
        :flux = tmean(flux_by_mlat)
        :n_time = ntime(flux_by_mlat)
        :nnan_count = count(!isnan, :flux)
    end
end

function get_flux_by_mlat(flux, mlat, timerange)
    flux_subset = tview(flux, timerange)
    mlat_subset = tview(mlat, timerange)
    return get_flux_by_mlat(flux_subset, mlat_subset)
end

end
