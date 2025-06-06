function dmsp_download(t::DateTime, id = 16, format = "hdf5")
    # <!-- '/opt/openmadrigal/madroot/experiments3/2020/dms/31dec20/dms_20201231_16e.001.hdf5' -->
    prefix = "/opt/openmadrigal/madroot/experiments3/"
    _dayofmonth = lpad(dayofmonth(t), 2, '0')
    _monthabbr = Dates.monthabbr(t) |> lowercase
    _year = string(Dates.year(t))
    _yearabbr = _year[3:4]
    _date = Dates.format(t, "yyyymmdd")
    filename = prefix * "$(_year)/dms/$(_dayofmonth)$(_monthabbr)$(_yearabbr)/dms_$(_date)_$(id)e.001.$format"
    return MadrigalWeb.download_file(filename)
end

function dmsp_download(timerange, id)
    files = map(timerange) do t
        dmsp_download(t, id)
    end
    return unique(files)
end

function dmsp_load(timerange, id, params::Vector{String})
    files = dmsp_download(timerange, id)
    return mapreduce(vcat, files) do file
        read2dimstack(file, params)
    end
end

function dmsp_load(timerange, id, param::String)
    files = dmsp_download(timerange, id)
    return if length(files) == 1
        read2dimarray(only(files), param, timerange)
    else
        vcat(read2dimarray.(files, param, Ref(timerange))...)
    end
end
