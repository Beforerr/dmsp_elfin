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
