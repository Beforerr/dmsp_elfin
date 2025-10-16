using Beforerr: fhist

const 𝐦 = (;
    mlat = :mlat => abs => "|MLAT| (degrees)",
    maxAE = :maxAE => "Max AE (nT)",
    E_c = :E_c => log10 => "Log E_c (keV)",
)

function _binedges(x)
    x == :mlat && return (50:0.5:80)
    x == :maxAE && return (0:50:1300)
    error("Unknown x: $x")
end

function _limits(x)
    x == :mlat && return (50, 80)
    x == :maxAE && return (0, 1300)
    error("Unknown x: $x")
end

_variable(x) = x isa Pair ? x[1] : x

fmt(from, to, i; leftclosed, rightclosed) = string(from + (to - from) / 2)
_mlat(s) = parse(Float64, String(s))
mlt_fmt(from, to, i; leftclosed, rightclosed) = string(Int(from + (to - from) / 2))

# Plot each variable for each AE range
# Create figure with 3 rows (one per AE range) and 6 columns (one per variable)
function plot_params_variation(f, df, vars; colorranges = (;), facet = :row)
    # Reorder MLT levels: 12→24, then 0→12
    levels!(df.mlt_bin, mlt_levels)
    # xtickformat = x -> string.(mod.(x .+ 12, 24))
    axis = (; xlabel = "MLT", ylabel = "MLAT")
    ae_bins = unique(df.maxAE_bin)

    ticks = WilkinsonTicks(4; k_max = 4)
    for (col_idx, varp) in enumerate(vars)
        var = _variable(varp)
        label = varp isa Pair ? varp[2] : string(varp)
        colorrange = get(colorranges, var, quantile(df[!, var], [0.02, 0.98]))
        for (row_idx, ae_bin) in enumerate(ae_bins)
            tdf_subset = @rsubset(df, :maxAE_bin == ae_bin; view = true)
            values = tdf_subset[!, var]
            fp = facet == :row ? (row_idx, col_idx) : (col_idx, row_idx)
            ax = Axis(f[fp...]; axis...)
            # Create heatmap
            x = tdf_subset.mlt_bin
            y = tdf_subset.mlat_bin
            # mlt_indices = mod.(x .+ 12, 24)
            mlt_indices = levelcode.(x)
            mlat_indices = _mlat.(y)
            hm = heatmap!(ax, mlt_indices, mlat_indices, values; colorrange)

            # Customize x-axis ticks for MLT
            ax.xticks = (1:length(mlt_levels), string.(CategoricalArrays.levels(x)))
            # Add colorbar
            cb_pos = facet == :row ? (row_idx - 1, col_idx) : (col_idx, row_idx - 1)
            vertical = !(facet == :row)
            row_idx == 1 && Colorbar(f[cb_pos...], hm; label, vertical = false, ticks)
            # Add AE range label on the left
            ae_pos = facet == :row ? (row_idx, 0) : (0, row_idx)
            col_idx == 1 && Label(f[ae_pos...], "AE: $(ae_bin)", rotation = π / 2, tellheight = false, font = :bold)

            row_idx != length(ae_bins) && hidexdecorations!(ax)
            col_idx != 1 && hideydecorations!(ax)
        end
        # plt = data(tdf_spatial) * mapping(:mlt_bin, :mlat_bin, var) * visual(Heatmap) * mapping(row=:maxAE_bin)
    end
    return
end

function plot_empirical_conditional_density(f, df, vars; mlt_idxs = 1:2:12)
    ssdf = @chain df begin
        @groupby(:mlt_bin, :mlat_bin)
        combine(vars .=> sum; renamecols = false)
    end

    foreach(enumerate(vars)) do (i, var)
        # emp = empirical_conditional_density(sdf, var)
        xlabel = string(var)
        hists = ssdf[!, var]
        limits = (0, quantile(_maximum.(hists), 0.95))
        idx = 1
        for (key, subdf) in pairs(@groupby(ssdf, :mlt_bin)[mlt_idxs])
            tdf = @rsubset(subdf, nentries($var) >= 8)
            matrix_data = extract_matrix(tdf, var)
            ax = Axis(f[idx, i]; xlabel, ylabel = "|MLAT| (deg)")
            heatmap!(ax, matrix_data; colorrange = limits)
            hists = tdf[!, var]
            mlats = _mlat.(tdf.mlat_bin)
            means = mean.(hists)
            stds = std.(hists)
            lines!(ax, means, mlats, color = :white, linewidth = 2)
            lines!(ax, means .+ stds, mlats, color = :white, linewidth = 1.5, linestyle = :dash)
            lines!(ax, means .- stds, mlats, color = :white, linewidth = 1.5, linestyle = :dash)
            ylims!(ax, (49, 81))
            idx != length(mlt_idxs) && hidexdecorations!(ax)
            i != 1 && hideydecorations!(ax)
            i == length(vars) && Label(f[idx, i + 1], "MLT = $(key.mlt_bin)", rotation = π / 2, tellheight = false)
            idx += 1
        end
        Colorbar(f[0, i]; limits, label = "n", vertical = false)
    end
    return f
end

function plot_parameter_distributions_aog(df, xsym = :mlat; bins = 60, normalization = :column, x_binedges = nothing)
    plot_df = @chain df begin
        @rtransform @astable begin
            :where = :mlat > 0 ? "North" : "South"
            :E_c = E_c(m1_len(:model))
            :γ = γ_len(m1_len(:model))
            :κ = κ_len(m2_len(:model))
            :E_c2 = E_c(m2_len(:model))
            :log10_E_c2 = log10(:E_c2)
        end
        @rsubset(isfinite(:E_c) & isfinite(:γ) && :E_c < 100, -10 < :κ < 15)
    end

    facet = mapping(layout = :where)
    base = data(plot_df) * facet

    x = getproperty(𝐦, xsym)
    vis = visual(colormap = :plasma, alpha = 0.8)

    # MLAT vs E_c density plot
    x_binedges = @something x_binedges _binedges(xsym)
    x_limits = extrema(x_binedges)

    plot_density = normalization != :column

    ec_plot = base * mapping(x, 𝐦.E_c)
    ec_plot *= if plot_density
        AoG.density(npoints = bins, datalimits = (x_limits, (-4, 2))) * vis
    else
        fhist(; npoints = bins, binedges = (x_binedges, -4:0.2:2), normalization) * vis
    end

    # MLAT vs γ density plot
    gamma_plot = base * mapping(x, :γ => "γ")
    gamma_plot *= if plot_density
        AoG.density(npoints = bins, datalimits = (x_limits, (5, -10))) * vis
    else
        fhist(; npoints = bins, binedges = (x_binedges, -10:0.5:5), normalization) * vis
    end

    # MLAT vs κ density plot
    κ_plot = base * mapping(x, :κ => "κ")
    κ_plot *= if plot_density
        AoG.density(npoints = bins, datalimits = (x_limits, (0, 12))) * vis
    else
        fhist(; npoints = bins, binedges = (x_binedges, 0:0.5:12), normalization) * vis
    end

    # MLAT vs E_c2 density plot
    ec2_plot = base * mapping(x, :log10_E_c2 => "Log E_c2 (keV)")
    ec2_plot *= if plot_density
        AoG.density(npoints = bins, datalimits = (x_limits, (-1, 2))) * vis
    else
        fhist(; npoints = bins, binedges = (x_binedges, -1:0.1:2), normalization) * vis
    end

    # Create figure with subplots
    f = Figure(size = (1200, 800))

    # Draw plots
    draw!(f[1, 1], ec_plot; facet = (; linkxaxes = :none))
    draw!(f[2, 1], gamma_plot; facet = (; linkxaxes = :none))
    draw!(f[3, 1], κ_plot; facet = (; linkxaxes = :none))
    draw!(f[4, 1], ec2_plot; facet = (; linkxaxes = :none))

    return f
end


function plot_all_means_by_ae_bin(df; binedges = 0:20:1000)
    # Extract all parameter values
    E_c_vals = E_c.(m1.(df.model))
    γ_vals = γ_len.(m1.(df.model))
    κ_vals = κ_len.(m2.(df.model))
    E_c2_vals = E_c.(m2.(df.model))
    mlat_vals = df.mlat
    Δmlt_vals = df.Δmlt

    # Create bins
    bin_indices = searchsortedlast.(Ref(binedges), df.maxAE)

    # Helper function to calculate stats for a parameter
    function calc_bin_stats(vals, pred = isfinite)
        means, stds, bin_centers, counts = Float64[], Float64[], Float64[], Int[]
        for i in 1:(length(binedges) - 1)
            mask = bin_indices .== i
            bin_vals = filter(pred, vals[mask])
            if !isempty(bin_vals)
                push!(means, mean(bin_vals))
                push!(stds, std(bin_vals))
                push!(bin_centers, (binedges[i] + binedges[i + 1]) / 2)
                push!(counts, length(bin_vals))
            end
        end
        return bin_centers, means, stds, counts
    end

    # Calculate stats for all parameters
    bc_Ec, means_Ec, stds_Ec, counts_Ec = calc_bin_stats(E_c_vals)
    bc_γ, means_γ, stds_γ, counts_γ = calc_bin_stats(γ_vals)
    bc_κ, means_κ, stds_κ, counts_κ = calc_bin_stats(κ_vals)
    bc_Ec2, means_Ec2, stds_Ec2, counts_Ec2 = calc_bin_stats(E_c2_vals)

    # Create figure with multiple panels
    f = Figure(size = (1400, 1000))

    # E_c panel
    ax1 = Axis(f[1, 1], xlabel = "AE (nT)", ylabel = "Mean E_c (keV)")
    errorbars!(ax1, bc_Ec, means_Ec, stds_Ec, whiskerwidth = 8)
    scatter!(ax1, bc_Ec, means_Ec, markersize = 12, color = :steelblue)
    ylims!(ax1, 0, 10)

    # γ panel
    ax2 = Axis(f[1, 2], xlabel = "AE (nT)", ylabel = "Mean γ")
    errorbars!(ax2, bc_γ, means_γ, stds_γ, whiskerwidth = 8)
    scatter!(ax2, bc_γ, means_γ, markersize = 12, color = :coral)
    ylims!(ax2, -10, 5)

    # κ panel
    ax3 = Axis(f[2, 1], xlabel = "AE (nT)", ylabel = "Mean κ")
    errorbars!(ax3, bc_κ, means_κ, stds_κ, whiskerwidth = 8)
    scatter!(ax3, bc_κ, means_κ, markersize = 12, color = :seagreen)
    # ylims!(ax3, 0, 12)

    # E_c2 panel
    ax4 = Axis(f[2, 2], xlabel = "AE (nT)", ylabel = "Mean E_c2 (keV)")
    errorbars!(ax4, bc_Ec2, means_Ec2, stds_Ec2, whiskerwidth = 8)
    scatter!(ax4, bc_Ec2, means_Ec2, markersize = 12, color = :purple)
    ylims!(ax4, 0, 30)

    # Count panel
    ax6 = Axis(f[3, 1:2], xlabel = "AE (nT)", ylabel = "Number of observations")
    barplot!(ax6, bc_Ec, counts_Ec, color = :gray, strokewidth = 1, strokecolor = :black)
    return f
end

"""
    plot_superposed_spectra!(fig_grid, df; mlat_range, mlt_range, n_sample=100, E_range=(0.03, 1000), n_points=200)

Plot superposed energy spectra for data within specified MLat and MLT ranges, separated by AE bins.

Creates a multi-panel visualization showing:
- Individual model spectra (semi-transparent lines)
- Median spectrum across all samples (red line)
- 25-75 percentile range (shaded band)

# Arguments
- `fig_grid`: Figure grid position (e.g., `f[1, :]` for row 1 of figure `f`)
- `df`: DataFrame containing fitted models with columns `:model`, `:mlat`, `:mlt_elx`, `:maxAE_bin`

# Keyword Arguments
- `mlat_range::Tuple{Real,Real}`: MLat range to filter (degrees)
- `mlt_range::Tuple{Real,Real}`: MLT range to filter (hours)
- `n_sample::Int=100`: Maximum number of spectra to plot per AE bin
- `show_legend::Bool=true`: Show legend on first panel
- `add_label::Bool=true`: Add label to the figure

# Dependencies
Requires StatsBase (for `sample`, `median`, `quantile`) and DataFramesMeta (for `@chain`, `@rsubset`)
to be available in the calling environment.

# Examples
```julia
# Single row with three AE bins
f = Figure(size=(1500, 500))
plot_superposed_spectra!(f[1, :], sdf; mlat_range=(70, 72), mlt_range=(5, 7))
```
"""
function plot_superposed_spectra!(
        fig_grid, df, E_grid;
        mlat_range,
        mlt_range,
        n_sample = 100,
        show_legend = true,
        add_label = true
    )
    fig_grid = GridLayout(fig_grid)
    # Filter data within specified MLat and MLT ranges
    filtered_df = @rsubset(
        df,
        mlat_range[1] <= abs(:mlat) <= mlat_range[2],
        mlt_range[1] <= :mlt_elx <= mlt_range[2]
    )

    # Get unique AE bins
    ae_bins = filter(!ismissing, sort(unique(filtered_df.maxAE_bin)))
    axs = []

    for (i, ae_bin) in enumerate(ae_bins)
        # Filter by AE bin
        ae_filtered = @rsubset(filtered_df, :maxAE_bin == ae_bin)
        n_total = nrow(ae_filtered)

        # Sample up to n_sample fittings
        n_plot = min(n_sample, n_total)
        if n_plot > 0
            sample_indices = sample(1:n_total, n_plot; replace = false)
            sample_df = ae_filtered[sample_indices, :]

            # Compute kappa statistics
            κ_values = sample_df.κ
            κ_median = round(median(κ_values), digits = 1)
            κ_q25 = round(quantile(κ_values, 0.25), digits = 1)
            κ_q75 = round(quantile(κ_values, 0.75), digits = 1)

            # Create axis with comprehensive title using LaTeX
            title_text = L"AE: %$(ae_bin),\, n: %$n_total,\, κ: %$(κ_median)_{%$(κ_q25)}^{%$(κ_q75)}"

            ax = Axis(
                fig_grid[1, i]; xlabel = "Energy (keV)", ylabel = YLabel.nflux,
                xscale = log10, yscale = log10,
                title = title_text
            )

            # Plot each model's spectrum
            for row in eachrow(sample_df)
                model = row.model
                flux_values = model.(E_grid)
                lines!(ax, E_grid, flux_values; alpha = 0.3, color = (:steelblue, 0.3))
            end

            # Compute and plot median spectrum
            flux_matrix = zeros(length(E_grid), nrow(sample_df))
            for (j, row) in enumerate(eachrow(sample_df))
                flux_matrix[:, j] .= row.model.(E_grid)
            end
            median_flux = vec(median(flux_matrix, dims = 2))
            lines!(ax, E_grid, median_flux; color = :red, linewidth = 3, label = "Median")

            # Compute and plot percentiles
            p25_flux = quantile.(eachrow(flux_matrix), 0.25)
            p75_flux = quantile.(eachrow(flux_matrix), 0.75)
            band!(ax, E_grid, p25_flux, p75_flux; color = (:red, 0.2), label = "25-75%")
            ylims!(ax, 1.0e1, 1.0e12)

            # Hide y decorations for columns after the first
            i != 1 && hideydecorations!(ax; ticks = false)

            (i == 1 && show_legend) && axislegend(ax; position = :rt)

            # Compute flux statistics (kappa component: model2)
            J_κ = median(sample_df.J2)  # Number flux (kappa)
            JE_κ = median(sample_df.JE2)  # Energy flux (kappa)

            # Compute total flux statistics (both components)
            J_total = median(sample_df.J1 .+ sample_df.J2)
            JE_total = median(sample_df.JE1 .+ sample_df.JE2)
            text = L"%$(latexxx((; J_κ, JE_κ)))\\%$(latexxx((; J_total, JE_total)))"
            text!(ax, 0.05, 0.05; text, space = :relative)

            push!(axs, ax)
        end
    end

    add_label && Label(
        fig_grid[0, 1:end], "Energy Spectra: MLat ∈ $(mlat_range), MLT ∈ $(mlt_range)",
        fontsize = 16, font = :bold
    )

    rowgap!(fig_grid, 2)

    return axs
end
