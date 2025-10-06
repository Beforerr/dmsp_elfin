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
import CDFDatasets as CDF

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
    return lazy"$BASE_URL/$probe/$level/$instrument/fast/$path1/{Y}/$(probe)_$(level)_$(datatype)_{Y}{M:02d}{D:02d}_$version.cdf"
end


function load(trange, probe = "a"; update::Bool = false, kw...)
    trange = DateTime.(trange)
    url = apply_date_format(build_url(; probe, kw...), trange[1])
    file = RemoteFile(url; dir = "elfin_data")
    (!isfile(file.path) || update) && download(file.uri, file.path)
    return CDFDataset(file.path)
end

end
