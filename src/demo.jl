using AlgebraOfGraphics
using CDAWeb
const COLOR_RANGE = (8.0e2, 1.2e11)

include("plot.jl")

set_Z_theme!() = begin
    set_aog_theme!()
    update_theme!(;
        figure_padding = 2,
    )
end


function _standardize(da)
    return set(da, Ti => DateTime.(da.dims[1].val))
end

# temporary fix for heatmap
function _heatmap!(ax, data; kw...)
    x = data.dims[1].val
    y = data.dims[2].val
    z = replace(data.data, 0 => NaN)
    return heatmap!(ax, x, y, z; kw...)
end

function simplify_axis_datetime_ticks!(ax)
    orig_ticks = ax.xticks[]
    tick_texts = map(enumerate(orig_ticks[2])) do (i, t)
        i == 1 ? t : split(t, "T")[2] # replace text before "T" with space
    end
    return ax.xticks[] = orig_ticks[1], tick_texts
end

function quicklook(f, trange, elx_flux, elx_mlt, elx_mlat, dmsp_flux, dmsp_mlt, dmsp_mlat; mlats = nothing, colorrange = (1.0e2, 2.0e11), colormap = :turbo, add_ae = true)
    tvars = map((elx_flux.prec, dmsp_flux)) do f
        tview(f, trange)
    end

    non_flux_axs = let layout = GridLayout(f[1, 1])

        ae = DimArray(CDAWeb.get_data("OMNI_HRO_1MIN/AE_INDEX", trange...)) |> _standardize
        ae_ax = Axis(layout[1, 1]; ylabel = "AE (nT)")
        scatterlines!(ae_ax, ae; label = "AE")

        markersize = 4
        label_kw = (; markersize = 8)
        mlt_ax = Axis(layout[2, 1]; ylabel = "MLT")
        scatter!(mlt_ax, _standardize(elx_mlt); label = "ELFIN" => label_kw, markersize)
        scatter!(mlt_ax, _standardize(dmsp_mlt); label = "DMSP" => label_kw, markersize)

        mlat_ax = Axis(layout[3, 1]; ylabel = "MLAT")
        scatter!(mlat_ax, _standardize(elx_mlat); label = "ELFIN" => label_kw, markersize)
        scatter!(mlat_ax, _standardize(dmsp_mlat); label = "DMSP" => label_kw, markersize)

        Legend(layout[2:3, 1, Right()], mlat_ax; markersize = 15)
        rowgap!(layout, 2)
        (ae_ax, mlt_ax, mlat_ax)
    end

    hidexdecorations!.(non_flux_axs; grid = false)

    flux_layout = GridLayout(f[2, 1])
    axs = map(enumerate(tvars)) do (i, f)
        ax = Axis(flux_layout[i, 1]; yscale = log10, ylabel = 𝒀.E)
        _heatmap!(ax, f; colorscale = log10, colorrange, colormap)
        ax
    end
    hidexdecorations!.(axs[1:(end - 1)]; grid = false)

    Colorbar(flux_layout[1:2, 2]; limits = colorrange, scale = log10, label = 𝒀.nflux, colormap)
    rowgap!(flux_layout, 0)
    colgap!(flux_layout, 0)

    linkxaxes!(axs..., non_flux_axs...)
    rowgap!(f.layout, 2)
    xlims!(current_axis(), trange...)
    return f
end

function demo_plot(trange, df, elx_flux, elx_mlt, elx_mlat, dmsp_flux, dmsp_mlt, dmsp_mlat; mlats = nothing, colormap = :turbo, colorrange = (1.0e2, 2.0e11), add_ratios = true)
    mlat_limits = extrema(df.mlat)

    sdf = @rsubset(df, :Δmlt < 1)
    mlats = @something mlats rand(sdf.mlat, 2)
    _sdf = @rsubset(sdf, :mlat ∈ mlats)

    elx_flux_ratio = ELFIN.flux_ratio(elx_flux)

    f = Figure(; size = (1000, 800))

    let layout = GridLayout(f[1, 1])

        p1 = plot_flux_by_mlat(layout[1, 1], replace(tview(elx_flux.prec, trange), 0 => NaN), elx_mlat; colorrange, colormap)
        p2 = plot_flux_by_mlat(layout[2, 1], replace(dmsp_flux, 0 => NaN), dmsp_mlat; colorrange, colormap)
        Colorbar(layout[1:2, 2], p1.plot; label = 𝒀.nflux)

        p3 = plot_flux_by_mlat(layout[3, 1], tview(elx_flux_ratio, trange), elx_mlat)
        Colorbar(layout[3, 2], p3.plot; label = L"j_{prec}/j_{trap}")

        ax = Axis(layout[4, 1]; ylabel = "MLT")
        scatter!(ax, set_mlat_dim(elx_mlt, elx_mlat); label = "ELFIN")
        scatter!(ax, set_mlat_dim(dmsp_mlt, dmsp_mlat); label = "DMSP")
        Legend(layout[4, 2], ax; padding = (30, 0, 0, 0), tellheight = false, tellwidth = false)

        mlat_axes = (p1.axis, p2.axis, p3.axis, ax)
        vlines!.(mlat_axes, (mlats,); linestyle = :dash)
        xlims!.(mlat_axes, mlat_limits...)
        hidexdecorations!.(mlat_axes[1:3]; grid = false)
        rowgap!(layout, 1, 0)
        colgap!(layout, 0)
        rowsize!(layout, 4, Auto(0.6))
        linkxaxes!(mlat_axes...)
    end

    let layout = GridLayout(f[2, 1])

        axs = plot_spectra(layout[1, 1:2], _sdf)

        hlines!.(axs, FLUX_THRESHOLD; color = :black, linestyle = :dash)

        if add_ratios
            hidexdecorations!.(axs; grid = false, ticklabels = false)

            fg = layout[2, 1:2]
            whistler_threshold = 0.25

            axs = map(enumerate(eachrow(_sdf))) do (i, row)
                ax = Axis(fg[1, i]; xlabel = 𝒀.E, ylabel = L"j_{prec}/j_{trap}", xscale = log10, yscale = log10)
                scatterlines!.(ax, eachrow(tview(elx_flux_ratio, row.trange_elx)))
                # add text with mlt and mechanism
                text = "$(row.mechanism) $(row.spins_mechanism)"
                text!(ax, 0.02, 1; text, space = :relative, align = (:left, :top))
                xlims!(ax, 30, 350)
                ylims!(ax, 0.08, 10)
                i == 1 || hideydecorations!(ax; grid = false)
                ax
            end
            hlines!.(axs, ([whistler_threshold, FLCS_THRESHOLD[]],); color = :black, linestyle = :dash)
            rowgap!(layout, 0)
            colgap!(layout, 0)
            rowsize!(layout, 1, Auto(1.5))
            linkaxes!(axs...)
        end
        rowgap!(f.layout, 1, 0)
    end
    rowsize!(f.layout, 1, Auto(4))
    return f
end
