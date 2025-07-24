using Madrigal
export dmsp_get_aacgm, dmsp_load, dmsp_download

function dmsp_download(t::DateTime, id=16; format="hdf5")
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

dmsp_download(t::AbstractString, args...) = dmsp_download(DateTime(t), args...)

function dmsp_download(timerange, id, args...)
    files = map(timerange) do t
        dmsp_download(t, id, args...)
    end
    return unique(files)
end

function dmsp_load(timerange, id, params)
    files = dmsp_download(timerange, id)
    sts = map(files) do file
        read2dimstack(file, params, timerange)
    end
    return cat(sts...; dims=Ti)
end

function dmsp_load(timerange, id, param::String)
    files = dmsp_download(timerange, id)
    return if length(files) == 1
        read2dimarray(only(files), param, timerange)
    else
        vcat(read2dimarray.(files, param, Ref(timerange))...)
    end
end

function dmsp_get_aacgm(timerange, id)
    dmsp_stack = dmsp_load(timerange, id, ("gdlat", "glon", "gdalt"))
    times = parent(dims(dmsp_stack, Ti))
    return geod2aacgm.(dmsp_stack.gdlat, dmsp_stack.glon, dmsp_stack.gdalt, times)
end