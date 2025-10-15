module DMSP

using Dates
using Dates: AbstractTime
using Madrigal
using DimensionalData
using GeoCotrans
using HDF5
using JLD2

export load

include("hdf5.jl")

function download(timerange, id, args...)
    kindat = 10200 + id
    files = download_files(:DMSP, kindat, timerange...)
    return unique!(sort!(files))
end

function load(timerange, id, params)
    files = download(timerange, id)
    sts = map(files) do file
        read2dimstack(file, params, timerange)
    end
    return cat(sts...; dims = Ti)
end

function load(timerange, id, param::String)
    files = download(timerange, id)
    return if isempty(files)
        @warn "DMSP data not available for ID $id, timerange $timerange"
        nothing
    elseif length(files) == 1
        read2dimarray(only(files), param, timerange)
    else
        vcat(read2dimarray.(files, param, (timerange,))...)
    end
end

geod(timerange, id) = load(timerange, id, ("gdlat", "glon", "gdalt"))

using DimensionalData: YDim, @dim
@dim Energy YDim "Energy"

# use keV as the basic unit for energy dimension
# set the flux unit to 1/cm²/s/sr/MeV
# 2021-01-27T13:41:05 -  2021-01-27T13:41:09 for DMSP 16 there are some values way too low like around 1e-15
function flux(timerange, id; kw...)
    f = load(timerange, id, "el_d_flux"; kw...)
    isnothing(f) && return nothing

    f = set(f, Y => Energy(dims(f, Y).val .* 1.0e-3))
    f[f .< 1.0e-8] .= 0
    f .*= 1.0e6
    f.metadata[:description] = "Diff electron num flux (1/cm^2/s/sr/MeV)"
    f.metadata[:ylabel] = "Energy (keV)"
    f.metadata[:yscale] = log10
    return f
end
end
