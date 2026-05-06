include(joinpath(@__DIR__, "data.jl"))

let vars = [:JE_tot => log10 => L"(a) $J_E$: Energy flux", :J_tot => log10 => L"(b) $J$: Number flux", :R_JE_e30 => L"(c) $J_E^{>30 \text{keV}} / J_E$", :R_J_e30 => L"(d) $J^{>30 \text{keV}} / J$"], axis = m2axis, func = median
    f = Figure(; size = (720, 540))
    tdf = _make_tdf(gdf, func, vars)
    plot_params_variation(f, tdf, vars; axis)
    colgap!(f.layout, 4)
    rowgap!(f.layout, 4)
    easy_save("e30_flux_ratio_mlt_mlat_$(func)")
end
