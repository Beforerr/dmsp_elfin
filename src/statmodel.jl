using StatsBase: fit
using StatsBase

struct HistFunc{T} <: Function
    edges::T
end

(f::HistFunc{T})(x) where {T} = Hist1D(x; binedges = f.edges)

function _hist(vals, bins)
    return StatsBase.fit(Histogram, vals, bins)
end

function extract_matrix(df, var::Symbol)
    if nrow(df) == 0
        return nothing
    end

    # Sort by mlat_bin to ensure consistent ordering
    sort!(df, :mlat_bin)

    # Get dimensions from first histogram
    first_hist = df[1, var]
    p_centers = bincenters(first_hist)
    n_p = length(p_centers)
    n_mlat = nrow(df)

    # Build matrix
    density_matrix = zeros(Float64, n_p, n_mlat)
    mlat_bins = parse.(Float64, string.(df.mlat_bin))

    for (j, row) in enumerate(eachrow(df))
        density_matrix[:, j] .= row[var].bincounts
    end

    return DimArray(density_matrix, (Dim{var}(p_centers), Dim{:mlat}(mlat_bins)))
end

function empirical_conditional_density(
        df,
        var::Symbol;
        n_var_bins = 31,
        n_mlat_bins = 19,
        n_mlt_bins = 12,
        mlat_range = (50, 80)
    )
    # Prepare data
    work_df = @chain df begin
        @select($var, :mlat, :elfin_mlt, :maxAE)
        @rtransform(:abs_mlat = abs(:mlat))
        @rsubset!(mlat_range[1] < :abs_mlat <= mlat_range[2])
    end

    var_vals = work_df[!, var]
    abs_mlat_vals = work_df.abs_mlat
    mlt_vals = work_df.elfin_mlt

    # Create bin edges for the variable
    p_min, p_max = extrema(var_vals)
    p_span = max(p_max - p_min, 1.0e-3)
    p_edges = range(p_min - 0.05 * p_span, p_max + 0.05 * p_span; length = n_var_bins)

    # Create bin edges for mlat
    mlat_min, mlat_max = extrema(abs_mlat_vals)
    mlat_span = max(mlat_max - mlat_min, 1.0e-3)
    mlat_edges = range(mlat_min - 0.05 * mlat_span, mlat_max + 0.05 * mlat_span; length = n_mlat_bins)

    # Create bin edges for mlt
    mlt_edges = range(0, 24; length = n_mlt_bins + 1)

    # Fit 3D histogram
    hist = fit(Histogram, (var_vals, abs_mlat_vals, mlt_vals), (p_edges, mlat_edges, mlt_edges))

    # Compute bin centers and widths
    p_centers, mlat_centers, mlt_centers = map(midpoints, hist.edges)
    widths_P = diff(hist.edges[1])

    # Compute conditional density p(var | mlat, mlt)
    cond_density = similar(hist.weights, Float64)

    @views for j in axes(cond_density, 2), k in axes(cond_density, 3)
        counts = hist.weights[:, j, k]
        total = sum(counts)
        if total > 0
            cond_density[:, j, k] .= counts ./ (total .* widths_P)
        else
            fill!(cond_density[:, j, k], 0.0)
        end
    end

    return (
        df = work_df,
        hist = hist,
        cond_density = cond_density,
        p_centers = p_centers,
        mlat_centers = mlat_centers,
        mlt_centers = mlt_centers,
        widths_P = widths_P,
        var = var,
    )
end
