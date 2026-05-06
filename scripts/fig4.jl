include(joinpath(@__DIR__, "data.jl"))

let axis = m2axis, func = median
    vars = [:E => L"(a) $\bar{E}$: Averaged energy", :κ => L"(b) $\kappa$: Kappa parameter"]
    colorranges = (; E = (0.1, 5.5), κ = (1.8, 6.1))
    f = Figure(; size = (450, 540))
    tdf = _make_tdf(gdf, func, vars)
    plot_params_variation(f, tdf, vars; axis, colorranges)
    colgap!(f.layout, 4)
    rowgap!(f.layout, 4)
    easy_save("key_params_mlt_mlat_$(func)")
end
