export get_elfin_precipitating_flux, get_elfin_gei

function standardize_pyspedas_dimarray(da)
    return set(da, Dim{:time} => Ti, Dim{:v_dim} => Y)
end

function get_elfin_precipitating_flux(trange, probe; level = "l2")
    epd_ds = elfin.epd(trange; probe, level)
    @info "availabel tvar names: $(keys(epd_ds))"
    elx_para = DimArray(getproperty(epd_ds, Symbol(:el, probe, :_pef_hs_nflux_para))) |> standardize_pyspedas_dimarray
    elx_anti = DimArray(getproperty(epd_ds, Symbol(:el, probe, :_pef_hs_nflux_anti))) |> standardize_pyspedas_dimarray
    elx_flux = DimArray(abs.(elx_para .- elx_anti))
    ds = DimStack((; para = elx_para, anti = elx_anti, prec = elx_flux))

    for da in layers(ds)
        replace!(da, get(da.metadata, "VALIDMIN", 0) => NaN)
        promote_cdf_attributes!(da.metadata)
        set_if_valid!(
            da.metadata,
            :yscale => log10, :ylabel => "Energy", :yunit => "keV",
            :scale => log10, :colorrange => (1.0e3, 1.0e8)
        )
    end

    return ds
end

function get_elfin_gei(trange, probe)
    elx_state = elfin.state(trange; probe)
    elx_gei = DimArray(getproperty(elx_state, Symbol(:el, probe, :_pos_gei)))
    return elx_gei = set(elx_gei, Dim{:time} => Ti)
end
