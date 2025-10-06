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
