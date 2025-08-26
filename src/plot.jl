using Makie: heatmap!, Axis
export plot_elfin_dmsp, plot_spectra
export plot_flux_by_mlat, plot_flux_by_mlat
using Makie
using Printf

import SpacePhysicsMakie
using SpacePhysicsMakie: set_if_valid!

baremodule YLabel
    nflux = "Flux (1/cm²/s/sr/MeV)"
    E = "Energy (keV)"
    B = "B (nT)"
    V = "V (km s⁻¹)"
    n = "n (cm⁻³)"
    ΔV = "ΔV (km/s)"
    T = "T (eV)"
    J = "J (nA m⁻²)"
end

const 𝒀 = YLabel


function set_dmsp_flux_opts!(ds)
    for f in (ds.el_d_flux, ds.ion_d_flux)
        flux_opts!(f)
    end
    return ds
end

set_mlat_dim(flux, mlat) = set(flux, Ti => X(mlat))
set_mlat_dim(flux, mlat, tspan) = set(tview(flux, tspan), Ti => X(tview(mlat, tspan)))

function plot_flux_by_mlat(f, flux; kw...)
    ax = Axis(f; yscale = log10, ylabel = 𝒀.E)
    plot = heatmap!(ax, flux; colorscale = log10, kw...)
    return Makie.AxisPlot(ax, plot)
end

function plot_flux_by_mlat(f, flux, mlat, args...; kw...)
    return plot_flux_by_mlat(f, set_mlat_dim(flux, mlat, args...); kw...)
end

function plot_elfin_dmsp(timerange, ids)
    dmsp_fluxs = map(ids) do id
        DMSP.load(timerange, id, "el_d_flux")
    end
    foreach(dmsp_fluxs) do flux
        set_if_valid!(
            flux.metadata,
            :yscale => log10,
            :scale => log10, :colorrange => (1.0e1, 1.0e6)
        )
        replace!(flux, 0 => NaN)
    end
    tvars = (elx_flux, elx_flux_perp, dmsp_fluxs..., elx_mlt, elx_aacgm.mlat, Δmlts, Δmlats)
    faxs = SpacePhysicsMakie.tplot(tvars, timerange...)
    ylims!.(faxs.axes[1:2], 70, 2.0e3)
    return faxs
end


"""
    plot_conjunction(event, dmsp_interp, elfin_interp, output_dir="plots")

Plot a conjunction event.
"""
function plot_conjunction(event, dmsp_interp, elfin_interp, output_dir = "plots")
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
    p = plot(layout = (2, 1), size = (800, 600), legend = :topright)

    # Plot MLT
    plot!(p[1], time, dmsp_interp.mlt[time_slice], label = "DMSP", color = :blue, lw = 2)
    plot!(p[1], time, elfin_interp.mlt[time_slice], label = "ELFIN", color = :red, lw = 2)

    # Add shaded area for conjunction
    vspan!(
        p[1], [Dates.value(event.start_time), Dates.value(event.end_time)],
        alpha = 0.2, color = :green, label = ""
    )

    ylabel!(p[1], "MLT (hours)")
    title!(p[1], "DMSP-ELFIN Conjunction Event")

    # Plot MLAT
    plot!(p[2], time, dmsp_interp.mlat[time_slice], label = "DMSP", color = :blue, lw = 2)
    plot!(p[2], time, elfin_interp.mlat[time_slice], label = "ELFIN", color = :red, lw = 2)

    # Add shaded area for conjunction
    vspan!(
        p[2], [Dates.value(event.start_time), Dates.value(event.end_time)],
        alpha = 0.2, color = :green, label = ""
    )

    ylabel!(p[2], "MLAT (degrees)")
    xlabel!(p[2], "Time (UTC)")

    # Add subtitle with event details
    start_str = Dates.format(event.start_time, "yyyy-mm-dd HH:MM:SS")
    duration_min = event.duration_seconds / 60
    title!(
        p, "Start: $start_str, Duration: $(round(duration_min, digits = 1)) minutes",
        subplot = 0
    )

    # Save figure
    filename = "conjunction_$(Dates.format(event.start_time, "yyyymmdd_HHMMSS")).png"
    filepath = joinpath(output_dir, filename)
    savefig(p, filepath)

    println("Saved plot to $filepath")
    return filepath
end

function plot_spectra!(ax, flux, flux_1, model)
    Emin = model.Emin
    scatterlines!(ax, flux)
    scatterlines!(ax, flux_1)

    energies = vcat(flux.dims[1].val, flux_1.dims[1].val)
    lines!(ax, energies, model.(energies); label = "Combined Model", linewidth = 3, color = :red)
    vlines!(ax, Emin; label = "Transition Energy", color = :black, linestyle = :dash)

    # Plot individual components with parameters in labels
    m1 = model.model1
    m2 = model.model2
    m1_label = "PLEC: A=$(@sprintf "%.0e" m1.A), γ=$(@sprintf "%.1e" m1.γ), Ec=$(@sprintf "%.0e" m1.E_c)keV"
    m2_label = string(m2)

    lines!(ax, energies, m1.(energies); label = m1_label, linestyle = :dot, color = :blue)
    E_high = energies[energies .>= Emin]
    lines!(ax, E_high, m2.(E_high); label = m2_label, linestyle = :dot, color = :green)
    return ax
end

function plot_spectra(f, args...; kw...)
    ax = Axis(f; xlabel = 𝒀.E, ylabel = 𝒀.nflux, xscale = log10, yscale = log10)
    plot_spectra!(ax, args...; kw...)
    return ax
end

function plot_spectra(f, df::DataFrame)
    axs = map(enumerate(eachrow(df))) do (i, row)
        ax = plot_spectra(f[1, i], row.flux, row.flux_1, row.model)
        i == 1 || hideydecorations!(ax; grid = false)
        ax
    end
    xlims!.(axs, 0.011, 20000)
    ylims!.(axs, 1.0e1, 1.0e12)
    return axs
end


export plot_example_fits, plot_parameters_variation
export plot_PowerLawExpCutoff_parameter_variation, plot_PowerLaw_parameter_variation, plot_SmoothBrokenPowerlaw_parameter_variation

# Plot : Example fitted spectra for selected MLATs
function plot_example_fits!(ax, df)
    successful_fits = filter(r -> r.success, df; view = true)
    # Plot a few example fits
    example_indices = [1, div(nrow(successful_fits), 2), nrow(successful_fits)]
    # color-coded by MLAT
    colors = [:blue, :red, :green]

    for (i, idx) in enumerate(example_indices)
        row = successful_fits[idx, :]  # Get corresponding row from dataframe

        # Get original data
        ff = vcat(row.flux.data, row.flux_1.data)
        ee = vcat(row.flux.dims[1].val, row.flux_1.dims[1].val)
        ff, ee = remove_nan(ff, ee)

        # Plot original data
        scatter!(
            ax, ee, ff, color = colors[i], alpha = 0.6,
            label = "MLAT $(round(row.mlat, digits = 1))°"
        )

        lines!(
            ax, ee, row.flux_modeled, color = colors[i], alpha = 0.6,
            label = "MLAT $(round(row.mlat, digits = 1))°"
        )
    end
    return
end


function plot_example_fits(f, args...)
    ax = Axis(
        f, xlabel = 𝒀.E, ylabel = 𝒀.nflux,
        xscale = log10, yscale = log10, title = "Fitted Spectra Examples"
    )
    plot_example_fits!(ax, args...)
    axislegend(ax, position = :lb)
    return ax
end

function plot_parameters_variation(f, mlats, models, n_points; scores = nothing)
    # Row 1: PowerLawExpCutoff parameters

    PowerLawExpCutoff_models = getindex.(models, 1)
    SmoothBrokenPowerlaw_models = getindex.(models, 2)
    Emins = getindex.(models, 3)

    plot_PowerLawExpCutoff_parameter_variation(f[1, 1][1:3, 1], mlats, PowerLawExpCutoff_models)
    plot_SmoothBrokenPowerlaw_parameter_variation(f[1, 2][1:5, 1], mlats, SmoothBrokenPowerlaw_models)
    plot_fit_parameters_variation(f[1, 1][4:5, 1], mlats, Emins, n_points; scores)
    return f
end

function plot_fit_parameters_variation(f, mlats, Emins, n_points; scores = nothing)
    xlabel = "MLAT"
    ax1 = Axis(f[1, 1]; xlabel, ylabel = "Energy Transition (keV)")
    scatterlines!(mlats, Emins, color = :gray)
    ax2 = Axis(f[2, 1]; xlabel, ylabel = "Number of Data Points")
    scatterlines!(mlats, n_points, color = :black)

    axs = [ax1, ax2]

    if !isnothing(scores)
        ax3 = Axis(f[3, 1]; xlabel, ylabel = "Fit Score")
        scatterlines!(mlats, scores, color = :green)
        push!(axs, ax3)
    end
    hidexdecorations!.(axs[1:(end - 1)]; grid = false)
    return f
end


function plot_PowerLawExpCutoff_parameter_variation(f, mlats, params)
    As = [p.A for p in params]
    γs = [p.γ for p in params]
    E_cs = [p.E_c for p in params]
    ax1 = Axis(f[1, 1]; ylabel = "Amplitude A", yscale = log10)
    scatterlines!(ax1, mlats, As, color = :blue)
    ax2 = Axis(f[2, 1]; ylabel = "Power Index γ")
    scatterlines!(ax2, mlats, γs, color = :red)
    ax3 = Axis(f[3, 1]; xlabel = "MLAT", ylabel = "Cutoff Energy E_c (keV)", yscale = log10)
    scatterlines!(ax3, mlats, E_cs, color = :green)
    hidexdecorations!.((ax1, ax2); grid = false)

    return f
end

function plot_PowerLaw_parameter_variation(f, mlats, params)
    As = [p.A for p in params]
    γs = [p.γ for p in params]
    ax1 = Axis(f[1, 1]; ylabel = "Amplitude A", yscale = log10)
    scatterlines!(ax1, mlats, As, color = :blue)
    ax2 = Axis(f[2, 1]; xlabel = "MLAT", ylabel = "Power Index γ")
    scatterlines!(ax2, mlats, γs, color = :red)
    hidexdecorations!(ax1; grid = false)

    return f
end

function plot_SmoothBrokenPowerlaw_parameter_variation(f, mlats, params)
    As = [p.A for p in params]
    γ1s = [p.γ1 for p in params]
    γ2s = [p.γ2 for p in params]
    Eb = [p.Eb for p in params]
    ms = [p.m for p in params]

    ax1 = Axis(f[1, 1]; ylabel = "Amplitude A", yscale = log10)
    scatterlines!(ax1, mlats, As, color = :purple)
    ax2 = Axis(f[2, 1]; ylabel = "Power Index γ₁")
    scatterlines!(ax2, mlats, γ1s, color = :orange)
    ax3 = Axis(f[3, 1]; ylabel = "Power Index γ₂")
    scatterlines!(ax3, mlats, γ2s, color = :cyan)
    ax4 = Axis(f[4, 1]; ylabel = "Break Energy Eb (keV)", yscale = log10)
    scatterlines!(ax4, mlats, Eb, color = :brown)
    ax5 = Axis(f[5, 1]; xlabel = "MLAT", ylabel = "Smoothness m")
    scatterlines!(ax5, mlats, ms, color = :teal)
    hidexdecorations!.((ax1, ax2, ax3, ax4); grid = false)
    return f
end
