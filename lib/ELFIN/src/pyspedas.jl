import PySPEDAS: promote_cdf_attributes!, get_data

function py_epd(trange, probe; level = "l2", collect = false, type = "nflux", fullspin = false, datatype = "pef", kw...)
    epd_ds = elfin.epd(string.(trange); probe, level, fullspin, type_ = string(type), datatype = string(datatype), kw...)
    if length(epd_ds) == 0
        return nothing
    end
    res = fullspin ? :fs : :hs
    @info "availabel tvar names: $(keys(epd_ds))"
    probe = "el$(probe)"
    base_var = join([probe, datatype, res, type], '_')
    omni = DimArray(PySPEDAS.get_data(Symbol(base_var, :_omni); collect)) |> standardize
    para = DimArray(PySPEDAS.get_data(Symbol(base_var, :_para); collect)) |> standardize
    anti = DimArray(PySPEDAS.get_data(Symbol(base_var, :_anti); collect)) |> standardize
    perp = DimArray(PySPEDAS.get_data(Symbol(base_var, :_perp); collect)) |> standardize
    ch0 = PySPEDAS.get_data(Symbol(base_var, :_ch0); collect).data
    ch1 = PySPEDAS.get_data(Symbol(base_var, :_ch1); collect).data
    ch2 = PySPEDAS.get_data(Symbol(base_var, :_ch2); collect).data
    ch3 = PySPEDAS.get_data(Symbol(base_var, :_ch3); collect).data
    return (; omni, para, anti, perp, ch0, ch1, ch2, ch3)
end
