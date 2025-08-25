module DMSP

using Dates
using Madrigal
using DimensionalData
using GeoCotrans
using HDF5
using JLD2

export load

include("hdf5.jl")

function download(t::DateTime, id = 16; format = "hdf5")
    # <!-- '/opt/openmadrigal/madroot/experiments3/2020/dms/31dec20/dms_20201231_16e.001.hdf5' -->
    prefix = "/opt/openmadrigal/madroot/experiments3/"
    _dayofmonth = lpad(dayofmonth(t), 2, '0')
    _monthabbr = Dates.monthabbr(t) |> lowercase
    _year = string(Dates.year(t))
    _yearabbr = _year[3:4]
    _date = Dates.format(t, "yyyymmdd")
    filename = prefix * "$(_year)/dms/$(_dayofmonth)$(_monthabbr)$(_yearabbr)/dms_$(_date)_$(id)e.001.$format"
    return Madrigal.download_file(filename)
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
    return if length(files) == 1
        read2dimarray(only(files), param, timerange)
    else
        vcat(read2dimarray.(files, param, Ref(timerange))...)
    end
end

geod(timerange, id) = load(timerange, id, ("gdlat", "glon", "gdalt"))

end
