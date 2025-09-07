module ELFIN

import PySPEDAS
using PySPEDAS.Projects: elfin
using DimensionalData
using PySPEDAS: promote_cdf_attributes!

export precipitating_flux, gei
using Dates

function standardize(x)
    return set(x, Dim{:time} => Ti, Dim{:v_dim} => Y)
end

# About 1.4s time resolution
function precipitating_flux(trange, probe; level = "l2")
    epd_ds = elfin.epd(string.(trange); probe, level)
    @info "availabel tvar names: $(keys(epd_ds))"
    elx_para = DimArray(epd_ds[Symbol(:el, probe, :_pef_hs_nflux_para)]) |> standardize
    elx_anti = DimArray(epd_ds[Symbol(:el, probe, :_pef_hs_nflux_anti)]) |> standardize
    
    elx_flux = abs.(elx_para .- elx_anti)
    ds = DimStack((; para = elx_para, anti = elx_anti, prec = elx_flux))

    for da in layers(ds)
        replace!(da, get(da.metadata, "VALIDMIN", 0) => NaN)
        promote_cdf_attributes!(da.metadata)
        # disable it since inconsistent with the size
        # PySPEDAS.resolve_metadata_dependencies!(da)
    end
    return ds
end

function gei(trange, probe)
    elx_state = elfin.state(string.(trange); probe)
    name = Symbol(:el, probe, :_pos_gei)
    return DimArray(elx_state[name]) |> standardize
end

end