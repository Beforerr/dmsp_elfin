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

import PySPEDAS
using PySPEDAS.Projects: elfin
using DimensionalData
using Speasy.PythonCall
using Speasy: get_product
using Memoization
using Speasy: pycdfpp
import Speasy.speasy as spz
using CDFDatasets
import CDFDatasets as CDF

export precipitating_flux, gei, epd
using Dates

const BASE_URL = "https://data.elfin.ucla.edu"

include("utils.jl")
include("state.jl")
include("epd.jl")
include("pyspedas.jl")


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


function load(trange, probe = "a"; no_update::Bool = false, kw...)
    trange = DateTime.(trange)
    url = apply_date_format(build_url(; probe, kw...), trange[1])
    file = RemoteFile(url; dir = "elfin_data")
    if !isfile(file.path)
        download(file.uri, file.path)
    end
    if no_update
    end
    return CDFDataset(file.path)
end

# About 1.4s time resolution
function precipitating_flux(trange, probe; level = "l2", kwargs...)
    ds = py_epd(trange, probe; level, kwargs...)

    elx_flux = abs.(elx_para .- elx_anti)
    ds = DimStack((; omni = elx_omni, para = elx_para, anti = elx_anti, prec = elx_flux))

    for da in layers(ds)
        meta = da.metadata
        # replace!(da, get(meta, "VALIDMIN", 0) => NaN)
        "CDF" in keys(meta) && promote_cdf_attributes!(meta)
        # remove any metadata that is a Python object (avoid save and load errors)
        for (k, v) in pairs(meta)
            ispy(v) && delete!(meta, k)
        end

        # disable it since inconsistent with the size
        # PySPEDAS.resolve_metadata_dependencies!(da)
    end
    return ds
end

end
