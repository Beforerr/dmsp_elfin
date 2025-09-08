module ELFIN

import PySPEDAS
using PySPEDAS.Projects: elfin
using DimensionalData
using PySPEDAS: promote_cdf_attributes!, get_data
using PySPEDAS.PythonCall
using Speasy: get_product
using Memoization
using Printf
using Speasy: pycdfpp
import Speasy.speasy as spz

export precipitating_flux, gei
using Dates

function standardize(x)
    return set(x, Dim{:time} => Ti, Dim{:v_dim} => Y)
end

# About 1.4s time resolution
@memoize function precipitating_flux(trange, probe; level = "l2", collect = true)
    epd_ds = elfin.epd(string.(trange); probe, level)
    @info "availabel tvar names: $(keys(epd_ds))"
    elx_para = DimArray(get_data(Symbol(:el, probe, :_pef_hs_nflux_para); collect)) |> standardize
    elx_anti = DimArray(get_data(Symbol(:el, probe, :_pef_hs_nflux_anti); collect)) |> standardize

    elx_flux = abs.(elx_para .- elx_anti)
    ds = DimStack((; para = elx_para, anti = elx_anti, prec = elx_flux))

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

function gei(probe, t0, t1)
    probe = "el$(probe)"
    url_pattern = "https://data.elfin.ucla.edu/$probe/l1/state/defn/{Y}/$(probe)_l1_state_defn_{Y}{M:02d}{D:02d}_v\\d+.cdf"
    var = get_product(
        url_pattern, "$(probe)_pos_gei", t0, t1;
        use_file_list = true,
    )
    return DimArray(var)
end

function gei(trange, probe)
    return gei(probe, trange[1], trange[2])
end

function build_url(probe, level, instrument, datatype, fluxtype, res)
    return lazy"$probe/$level/$instrument/fast/electron/{Y}/$(probe)_$(level)_epdef_{Y}{M:02d}{D:02d}_v\\d+.cdf"
end

function apply_date_format(pattern, date)
    return replace(
        pattern,
        "{Y}" => year(date),
        "{M:02d}" => string(month(date); pad = 2),
        "{D:02d}" => string(day(date); pad = 2)
    )
end

function epd(trange, probe; level = "l2")
    probe = "el$(probe)"
    instrument = "epd"
    datatype = "e"
    fluxtype = "nflux"
    res = "hs"
    base_url = "https://data.elfin.ucla.edu"

    url_pattern = "$(base_url)/$(probe)/$level/$instrument/fast/electron/{Y}/$(probe)_$(level)_epdef_{Y}{M:02d}{D:02d}_v\\d+.cdf"

    file_path = apply_date_format(build_url(probe, level, instrument, datatype, fluxtype, res), trange[1])

    return pycdfpp.load(spz.core.any_files.any_loc_open("$base_url/$file_path").read())
    # name = @printf("p%sf_%s_Epat_%s", datatype, res, fluxtype)

    # _build_url
    # var = get_product(
    #     url_pattern, name, trange[1], trange[2];
    #     use_file_list = true,
    # )
    # return DimArray(var)
end


# Use pyspedas to get gei (much slower)
function pyspedas_gei(trange, probe)
    elx_state = elfin.state(string.(trange); probe)
    name = Symbol(:el, probe, :_pos_gei)
    return DimArray(elx_state[name]) |> standardize
end

end
