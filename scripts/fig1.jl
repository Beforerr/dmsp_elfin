include("setup.jl")

using ELFINData
using GeoCotrans, GeoAACGM
using Beforerr: add_labels!, add_topbar!, easy_save
import Beforerr

Beforerr.DEFAULT_FORMATS = [:png, :pdf]

include("../src/demo.jl")
set_Z_theme!()

# id, trange = 17, [DateTime("2021-12-01T22:19:00"), DateTime("2021-12-01T22:28")]
probe, id, trange = "a", 16, [DateTime("2022-07-02T05:00"), DateTime("2022-07-02T05:06")]

t0, t1 = trange
config = (; probe, id, trange, t0, t1)

Δt = Minute(10)

elx_flux = ELFINData.epd_spectral(trange; probe) |> permutedims
elx_gei_trange = extend(trange, Second(1))
elx_gei = tview(
    DimArray(probe == "a" ? ELA_POS_GEI(elx_gei_trange) : ELB_POS_GEI(elx_gei_trange)),
    elx_gei_trange...
)
elx_aacgm = gei2aacgm(elx_gei)
elx_mlat = _mlat(elx_aacgm)
elx_geo = gei2geo(elx_gei)
elx_mlt = GeoCotrans.get_mlt(elx_geo)

dmsp_flux = DMSP.flux(extend(trange, Δt), id)
dmsp_mlt, dmsp_mlat = get_mlt_mlat(id, extend(trange, Δt))
zoom_tr = timerange(tview(elx_flux.para, trange...))

elx_mlt = rebuild(elx_mlt; metadata = Dict(:label => uppercase("ELFIN-$probe")))
dmsp_mlt = rebuild(dmsp_mlt; metadata = Dict(:label => uppercase("DMSP-F$id")))

tdf = workload(trange, (id,), elx_flux, elx_gei)

f = let mlats = [-70, -68.5]
    f = Figure(; size = (900, 700))
    flux1 = elx_flux.prec
    flux2 = dmsp_flux
    demo_plot(
        f, zoom_tr, tdf, flux1, elx_mlt, elx_mlat, flux2, dmsp_mlt, dmsp_mlat; mlats, colormap = :cividis, mlt_offset = -3, legend = nothing
    )
    axs = filter(x -> x isa Axis, f.content)
    ax1 = axs[1]
    xlims!(ax1, -71, -52)
    ylims!(ax1, 0.8, 6.8)
    add_labels!(; position = Right(), padding = (0, 0, 0, 0))
    intervals = [(-70.5, -67)]
    add_topbar!(ax1, intervals)
    vspan!.(axs[4:5], flux2.dims[2][end], flux1.dims[2][1]; color = (:gray, 0.15))  # 30-50 keV energy gap between DMSP and ELFIN
    # easy_save("flux_with_fit")
    # easy_save(savename("flux_with_fit", config))
    f
end
