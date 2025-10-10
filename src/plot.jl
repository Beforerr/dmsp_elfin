import SpacePhysicsMakie
using SpacePhysicsMakie: set_if_valid!
using DmspElfinConjunction: TwoStepModel
using SpectralModels: math_show
import SpectralModels as SM
using UnPack
using TimeseriesUtilities: times, tinterp
using DimensionalData: set

# Use CairoMakie in CI, GLMakie otherwise
using GLMakie

# Set better default colormap for wide dynamic range (10^3 to 10^11)
set_theme!(colormap = :turbo)  # Excellent perceptual uniformity and high contrast

baremodule YLabel
    using LaTeXStrings
    nflux = "Flux (1/cm²/s/sr/MeV)"
    flux_ratio = L"j_{prec}/j_{trap}"
    E = "Energy (keV)"
    A = "Amplitude A"
    γ = γ1 = γ2 = "Power Index γ"
    E_c = "Cutoff Energy E_c (keV)"
    Eb = "Break Energy Eb (keV)"
    κ = "Kappa κ"
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

const COLORRANGE = Ref((1.0e3, 1.0e11))

export COLORRANGE

function set_mlat_dim(A, mlat)
    A_times = times(A)
    x = if times(mlat) == A_times
        mlat.data
    else
        _mlat = tinterp(mlat, A_times)
        _mlat.data
    end
    return set(A, Ti => X(x))
end

function plot_flux_by_mlat(f, flux; kw...)
    ax = Axis(f; yscale = log10, ylabel = 𝒀.E)
    f_kw = Dict()
    if haskey(flux.metadata, :colorrange)
        f_kw[:colorrange] = flux.metadata[:colorrange]
    end
    plot = heatmap!(ax, flux; colorscale = log10, f_kw..., kw...)
    return Makie.AxisPlot(ax, plot)
end

function plot_x_by_mlat(f, flux, mlat, args...; kw...)
    return plot_flux_by_mlat(f, set_mlat_dim(flux, mlat, args...); kw...)
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
    plot_conjunction(event, dmsp_interp, elfin_interp)

Plot a conjunction event.
"""
function plot_conjunction(event, dmsp_interp, elfin_interp)
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
    return p
end

function plot_spectra!(ax, energies, model)
    return lines!(ax, energies, model; label = SM.math_show(model), linestyle = :dot)
end

function plot_spectra!(ax, energies, model::TwoStepModel)
    Emin = model.Emin
    lines!(ax, energies, model.(energies); label = "Combined Model", linewidth = 2, color = :red)
    vlines!(ax, Emin; label = "Transition Energy", color = :black, linestyle = :dash)
    # axislegend(ax)

    # Plot individual components with parameters in labels
    l1 = plot_spectra!(ax, energies, model.model1)
    l2 = plot_spectra!(ax, energies, model.model2)
    axislegend(ax, [l1, l2], [l1.label, l2.label])
    return ax
end

function plot_spectra!(ax, flux, flux_1, model; plot_model = true)
    scatterlines!(ax, flux)
    scatterlines!(ax, flux_1)
    plot_model && begin
        energies = vcat(flux.dims[1].val, flux_1.dims[1].val)
        plot_spectra!(ax, energies, model)
    end
    return ax
end

function plot_spectra(f, args...; kw...)
    ax = Axis(f; xlabel = 𝒀.E, ylabel = 𝒀.nflux, xscale = log10, yscale = log10)
    plot_spectra!(ax, args...; kw...)
    return ax
end

function plot_spectra(f, df::DataFrame; kw...)
    axs = map(enumerate(eachrow(df))) do (i, row)
        ax = plot_spectra(f[1, i], row.flux_dmsp, row.flux_elx, row.model; kw...)
        i == 1 || hideydecorations!(ax; grid = false)
        ax
    end
    xlims!.(axs, 0.011, 20000)
    ylims!.(axs, 1.0e1, 1.0e12)
    linkyaxes!(axs...)
    linkxaxes!(axs...)
    return axs
end


export plot_parameters_variation, plot_flux_analysis
export plot_PowerLawExpCutoff_parameter_variation, plot_PowerLaw_parameter_variation, plot_SmoothBrokenPowerlaw_parameter_variation, plot_model_parameters, plot_parameters_variation_generic

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

"""
    plot_parameters_variation(f, mlats, models::Vector{<:SpectralModel}; colors=nothing)

Generic function to plot parameters of any spectral model type.
"""
function plot_parameters_variation(f, T, mlats, models)
    if isempty(models)
        return f
    end
    # Get field names from the first model
    field_names = fieldnames(T)
    # Create plots for each parameter
    axes = []
    for (i, fn) in enumerate(field_names)
        values = getfield.(models, fn)
        # Determine if we need log scale (for amplitude-like parameters)
        use_log = fn in (:A, :E_c, :Eb) || any(v -> v > 1000, values)

        ylabel = isdefined(YLabel, fn) ? getproperty(YLabel, fn) : string(fn)
        xlabel = "MLAT"

        ax = Axis(f[i, 1]; xlabel, ylabel, yscale = use_log ? log10 : identity)
        scatterlines!(ax, mlats, values)
        push!(axes, ax)
    end

    # Hide x decorations for all but last axis
    hidexdecorations!.(axes[1:(end - 1)]; grid = false)
    return f
end

"""
    plot_parameters_variation(f, TwoStepModel, mlats, models, n_points; scores=nothing)

Generalized version of plot_parameters_variation that works with any TwoStepModel types.

# Arguments
- `f`: Figure for plotting
- `n_points`: Number of data points for each fit
- `scores`: Optional fit quality scores

# Example  
```julia
models = [TwoStepModel(PowerLaw(...), KappaDistribution(...), 50.0), ...]
plot_parameters_variation(fig, mlats, models, n_points)
```
"""
function plot_parameters_variation(f, T, mlats, models, n_points; scores = nothing)
    # Extract model components
    model1s = [m.model1 for m in models]
    model2s = [m.model2 for m in models]
    Emins = [m.Emin for m in models]

    # Get model types
    M1_type = spectral_type(model1s)
    M2_type = spectral_type(model2s)

    # Calculate grid layout
    n_params_1 = length(fieldnames(M1_type))
    n_params_2 = length(fieldnames(M2_type))
    fit_params_rows = scores === nothing ? 2 : 3

    # Plot model1 parameters
    plot_parameters_variation(f[1:n_params_1, 1], mlats, model1s)

    # Plot model2 parameters
    plot_parameters_variation(f[1:n_params_2, 2], mlats, model2s)

    # Plot fit parameters (Emin, n_points, scores)
    fit_start_row = n_params_1 + 1
    plot_fit_parameters_variation(
        f[fit_start_row:(fit_start_row + fit_params_rows - 1), 1],
        mlats, Emins, n_points; scores
    )
    return f
end

spectral_type(T) = T
spectral_type(T::Union) = T.a == Nothing ? T.b : T.a
spectral_type(models::AbstractArray) = spectral_type(eltype(models))
plot_parameters_variation(f, mlats, models, args...; kw...) = plot_parameters_variation(f, spectral_type(models), mlats, models, args...; kw...)

dmsp_default_energy() = [0.03, 0.044, 0.065, 0.095, 0.139, 0.204, 0.3, 0.44, 0.646, 0.949, 1.392, 2.04, 3.0, 4.4, 6.46, 9.45, 13.9, 20.4, 30.0]
elfin_default_energy() = [63.25, 97.98, 138.56, 183.3, 238.12, 305.2, 385.16, 520.48, 752.99, 1081.67, 1529.71, 2121.32, 2893.96, 3728.61, 4906.12, 6500.0]
default_energies() = vcat(dmsp_default_energy(), elfin_default_energy())
