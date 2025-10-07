"""
Electron Losses and Fields Investigation (ELFIN)

# Instruments
- Fluxgate Magnetometer (FGM)
- Energetic Particle Detector (EPD)
- Magneto Resistive Magnetometer-a (MRMa)
- Magneto Resistive Magnetometer-i (MRMi)
- State data (state)
- Engineering data (ENG)
"""
module ELFIN
using DimensionalData
using Speasy: get_product
using Memoization
using CDFDatasets
using CDFDatasets: CDFDataset
import CDFDatasets.CommonDataModel as CDM
import CDFDatasets as CDF
using Downloads: request

export precipitating_flux, gei, epd
using Dates

const BASE_URL = "https://data.elfin.ucla.edu"

include("utils.jl")
include("state.jl")
include("epd.jl")

function build_url(; probe = "a", level = "l2", instrument = "epd", datatype = "pef", version = "v01")
    # Build URL pattern for ELFIN EPD data files
    !startswith(probe, "el") && (probe = "el$(probe)")
    if datatype == "pef"
        datatype = "epdef"
        path1 = "electron"
    elseif datatype == "pif"
        datatype = "epdif"
        path1 = "ion"
    end
    return "$BASE_URL/$probe/$level/$instrument/fast/$path1/{Y}/$(probe)_$(level)_$(datatype)_{Y}{M:02d}{D:02d}_$version.cdf"
end


function _download(url, output)
    mkpath(dirname(output))
    response = request(url; output)
    return if response.status == 200
        output
    else
        rm(output)
        @info "Remote file not available" url = url message = response.message
        nothing
    end
end

_download(file::RemoteFile) = _download(file.uri.uri, file.path)

function epd_download(trange, probe = "a"; update::Bool = false, dir = "elfin_data", kw...)
    dt = Day(1)
    t0 = floor(DateTime(trange[1]), dt)
    t1 = ceil(DateTime(trange[2]), dt)
    tranges = t0:dt:t1
    pattern = build_url(; probe, kw...)
    outputs = map(tranges) do ti
        url = format_date(pattern, ti)
        file = RemoteFile(url; dir)
        output = file.path
        (!isfile(output) || update) ? _download(file) : output
    end
    return filter!(!isnothing, outputs)
end

end
