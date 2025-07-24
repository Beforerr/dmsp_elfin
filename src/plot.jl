using Makie: heatmap!
export plot_elfin_dmsp, plot_x_mlat_flux!
export set_flux_opts!

using SPEDAS.TPlot: set_if_valid!

function set_flux_opts!(flux, range)
    set_if_valid!(flux.metadata, :yscale => log10, :scale => log10, :colorrange => range)
    replace!(flux, 0 => NaN)
    return flux
end

function set_dmsp_flux_opts!(ds)
    for f in (ds.el_d_flux, ds.ion_d_flux)
        flux_opts!(f)
    end
    return ds
end

function plot_x_mlat_flux!(ax, mlat, flux)
    y = flux.dims[2].val
    z = flux.data
    attrs = SPEDAS.TPlot.heatmap_attributes(flux)
    heatmap!(ax, mlat, y, z; attrs...)
end

function plot_elfin_dmsp(timerange, ids)
    dmsp_fluxs = map(ids) do id
        dmsp_load(timerange, id, "el_d_flux")
    end
    foreach(dmsp_fluxs) do flux
        set_if_valid!(flux.metadata,
            :yscale => log10,
            :scale => log10, :colorrange => (1e1, 1e6)
        )
        replace!(flux, 0 => NaN)
    end
    tvars = (elx_flux, elx_flux_perp, dmsp_fluxs..., elx_mlt, elx_aacgm.mlat, Δmlts, Δmlats)
    faxs = SPEDAS.tplot(tvars, timerange...)
    ylims!.(faxs.axes[1:2], 70, 2e3)
    faxs
end


"""
    plot_conjunction(event, dmsp_interp, elfin_interp, output_dir="plots")

Plot a conjunction event.
"""
function plot_conjunction(event, dmsp_interp, elfin_interp, output_dir="plots")
    # Create output directory if it doesn't exist
    mkpath(output_dir)

    # Extract event data
    start_idx = event.start_idx
    end_idx = event.end_idx

    # Add some margin for context
    margin = Int(min(300, round((end_idx - start_idx) * 0.5)))  # 5 minutes or 50% of event duration
    plot_start = max(1, start_idx - margin)
    plot_end = min(nrow(dmsp_interp), end_idx + margin)

    # Extract time slice for plotting
    time_slice = plot_start:plot_end
    time = dmsp_interp.timestamp[time_slice]

    # Create plot
    p = plot(layout=(2, 1), size=(800, 600), legend=:topright)

    # Plot MLT
    plot!(p[1], time, dmsp_interp.mlt[time_slice], label="DMSP", color=:blue, lw=2)
    plot!(p[1], time, elfin_interp.mlt[time_slice], label="ELFIN", color=:red, lw=2)

    # Add shaded area for conjunction
    vspan!(p[1], [Dates.value(event.start_time), Dates.value(event.end_time)],
        alpha=0.2, color=:green, label="")

    ylabel!(p[1], "MLT (hours)")
    title!(p[1], "DMSP-ELFIN Conjunction Event")

    # Plot MLAT
    plot!(p[2], time, dmsp_interp.mlat[time_slice], label="DMSP", color=:blue, lw=2)
    plot!(p[2], time, elfin_interp.mlat[time_slice], label="ELFIN", color=:red, lw=2)

    # Add shaded area for conjunction
    vspan!(p[2], [Dates.value(event.start_time), Dates.value(event.end_time)],
        alpha=0.2, color=:green, label="")

    ylabel!(p[2], "MLAT (degrees)")
    xlabel!(p[2], "Time (UTC)")

    # Add subtitle with event details
    start_str = Dates.format(event.start_time, "yyyy-mm-dd HH:MM:SS")
    duration_min = event.duration_seconds / 60
    title!(p, "Start: $start_str, Duration: $(round(duration_min, digits=1)) minutes",
        subplot=0)

    # Save figure
    filename = "conjunction_$(Dates.format(event.start_time, "yyyymmdd_HHMMSS")).png"
    filepath = joinpath(output_dir, filename)
    savefig(p, filepath)

    println("Saved plot to $filepath")
    return filepath
end