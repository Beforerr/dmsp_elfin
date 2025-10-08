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

function download(t::AbstractTime, id = 16; format = "hdf5")
    # <!-- '/opt/openmadrigal/madroot/experiments3/2020/dms/31dec20/dms_20201231_16e.001.hdf5' -->
    prefix = "/opt/openmadrigal/madroot/experiments3/"
    _dayofmonth = lpad(dayofmonth(t), 2, '0')
    _monthabbr = Dates.monthabbr(t) |> lowercase
    _year = string(Dates.year(t))
    _yearabbr = _year[3:4]
    _date = Dates.format(t, "yyyymmdd")
    filename = prefix * "$(_year)/dms/$(_dayofmonth)$(_monthabbr)$(_yearabbr)/dms_$(_date)_$(id)e.001.$format"
    return Madrigal.download_file(filename; download = (; update_period = 8))
end

download(t::AbstractString, args...) = download(DateTime(t), args...)

function download(timerange, id, args...)
    files = map(timerange) do t
        download(t, id, args...)
    end
    return unique(files)
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

    # Handle case when download returns nothing or empty array
    if isnothing(files) || isempty(files) || any(isnothing, files)
        @warn "DMSP data not available for ID $id, timerange $timerange"
        return nothing
    end

    return if length(files) == 1
        read2dimarray(only(files), param, timerange)
    else
        vcat(read2dimarray.(files, param, Ref(timerange))...)
    end
end

geod(timerange, id) = load(timerange, id, ("gdlat", "glon", "gdalt"))

# use keV as the basic unit for energy dimension
# set the flux unit to 1/cm²/s/sr/MeV
# 2021-01-27T13:41:05 -  2021-01-27T13:41:09 for DMSP 16 there are some values way too low like around 1e-15
function flux(timerange, id; kw...)
    f = load(timerange, id, "el_d_flux"; kw...)
    isnothing(f) && return nothing

    f = set(f, Y => dims(f, Y).val .* 1.0e-3)
    f[f .< 1.0e-8] .= 0
    f .*= 1.0e6
    return f
end
end
