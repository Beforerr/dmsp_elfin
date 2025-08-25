module ELFIN

import PySPEDAS
using PySPEDAS.Projects: elfin
using DimensionalData
using SpacePhysicsMakie: set_if_valid!
using PySPEDAS: promote_cdf_attributes!

# About 1.4s time resolution
function precipitating_flux(trange, probe; level = "l2")
    epd_ds = elfin.epd(trange; probe, level)
    @info "availabel tvar names: $(keys(epd_ds))"
    elx_para = DimArray(getproperty(epd_ds, Symbol(:el, probe, :_pef_hs_nflux_para)))
    elx_anti = DimArray(getproperty(epd_ds, Symbol(:el, probe, :_pef_hs_nflux_anti)))
    elx_energies = [63.0, 98.0, 139.0, 183.0, 238.0, 305.0, 385.0, 520.0, 753.0, 1082.0, 1530.0, 2121.0, 2894.0, 3729.0, 4906.0, 6500.0]
    dims = (Ti(parent(elx_para.dims[1])), Y(elx_energies))
    elx_para = DimArray(elx_para, dims; metadata = elx_para.metadata)
    elx_anti = DimArray(elx_anti, dims; metadata = elx_anti.metadata)

    elx_flux = abs.(elx_para .- elx_anti)
    ds = DimStack((; para = elx_para, anti = elx_anti, prec = elx_flux))

    for da in layers(ds)
        replace!(da, get(da.metadata, "VALIDMIN", 0) => NaN)
        promote_cdf_attributes!(da.metadata)
        set_if_valid!(
            da.metadata,
            :yscale => log10, :ylabel => "Energy", :yunit => "keV",
            :colorscale => log10, :colorrange => (1.0e3, 1.0e8)
        )
        # disable it since inconsistent with the size
        # PySPEDAS.resolve_metadata_dependencies!(da)
    end
    return ds
end

function gei(trange, probe)
    elx_state = elfin.state(trange; probe)
    elx_gei = DimArray(getproperty(elx_state, Symbol(:el, probe, :_pos_gei)))
    return set(elx_gei, Dim{:time} => Ti)
end


end