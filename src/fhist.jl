# https://github.com/MakieOrg/AlgebraOfGraphics.jl/blob/master/src/transformations/density.jl#L18

module FHistAog
using FHist
using CairoMakie.Makie
using AlgebraOfGraphics
using AlgebraOfGraphics: automatic, transformation, set, nested_extrema_finite, valid_options
using AlgebraOfGraphics: ProcessedLayer
export fhist

Base.@kwdef struct FHistAnalysis{D, K, B}
    datalimits::D = automatic
    npoints::Int = 200
    kernel::K = automatic
    binedges::B = automatic
    normalization::Symbol = :none
end

defaultdatalimits(positional) = map(nested_extrema_finite, Tuple(positional))

applydatalimits(f::Function, d) = map(f, d)
applydatalimits(limits::Tuple{Real, Real}, d) = map(_ -> limits, d)
applydatalimits(limits::Tuple, _) = limits

function _hist(data::NTuple{N}; kwargs...) where {N}
    N == 1 && return Hist1D(data...; kwargs...)
    N == 2 && return Hist2D(data; kwargs...)
    N == 3 && return Hist3D(data; kwargs...)
    throw(ArgumentError("Only 1D, 2D and 3D histograms are supported"))
end

function _density(vs::Tuple; normalization = :none, kwargs...)
    h = _hist(vs; kwargs...)
    z = bincounts(h)

    # Apply normalization
    if normalization == :column
        # Normalize each column (x-value) to sum to 1
        if ndims(z) == 2
            # For 2D histograms, normalize along the second dimension (columns)
            col_sums = sum(z, dims = 2)
            # Avoid division by zero
            col_sums[col_sums .== 0] .= 1
            z = z ./ col_sums
        elseif ndims(z) == 1
            # For 1D histograms, normalize entire array
            total = sum(z)
            if total > 0
                z = z ./ total
            end
        else
            @warn "Column normalization not implemented for $(ndims(z))D histograms"
        end
    elseif normalization == :probability
        # Standard probability normalization (entire histogram sums to 1)
        total = sum(z)
        if total > 0
            z = z ./ total
        end
    elseif normalization == :density
        # Density normalization (integral equals 1)
        edges = binedges(h)
        if length(edges) == 2
            # 2D histogram
            dx = diff(edges[1])
            dy = diff(edges[2])
            area = dx * dy'  # Broadcasting to create area matrix
            total_area = sum(z .* area)
            if total_area > 0
                z = z ./ total_area
            end
        elseif length(edges) == 1
            # 1D histogram
            dx = diff(edges[1])
            total_area = sum(z .* dx)
            if total_area > 0
                z = z ./ total_area
            end
        end
    end
    # :none normalization returns raw counts

    return (binedges(h)..., z)
end

function (d::FHistAnalysis)(input::ProcessedLayer)
    datalimits = d.datalimits === automatic ? defaultdatalimits(input.positional) : d.datalimits
    options = valid_options(; d.binedges, d.normalization)
    output = map(input) do p, n
        @info n
        return _density(Tuple(p); pairs(n)..., pairs(options)...), (;)
    end
    N = length(input.positional)
    labels = set(input.labels, N + 1 => "pdf")
    plottypes = [LinesFill, Heatmap, Volume]
    default_plottype = plottypes[N]
    plottype = Makie.plottype(input.plottype, default_plottype)
    return ProcessedLayer(output; plottype, labels)
end

"""
    fhist(; datalimits=automatic, npoints=200, kernel=automatic, binedges=automatic, normalization=:none)

Create a histogram using FHist with various normalization options.

# Arguments
- `datalimits`: Range for histogram calculation (defaults to data extrema)
- `npoints`: Number of points for visualization
- `kernel`: Kernel for density estimation (if applicable)
- `binedges`: Custom bin edges (automatic by default)
- `normalization`: Normalization method:
  - `:none` - Raw counts (default)
  - `:column` - Each column (x-value) sums to 1 
  - `:probability` - Entire histogram sums to 1
  - `:density` - Integral equals 1 (proper density)

# Examples
```julia
# Raw counts histogram
data(df) * mapping(:x, :y) * fhist()

# Column-normalized (each x-value column sums to 1)
data(df) * mapping(:x, :y) * fhist(normalization=:column)

# Probability histogram (entire histogram sums to 1)
data(df) * mapping(:x, :y) * fhist(normalization=:probability)

# Density histogram (integral equals 1)
data(df) * mapping(:x, :y) * fhist(normalization=:density)
```

The `:column` normalization is particularly useful for conditional probability distributions,
where you want to see the probability distribution of y given each x value.
"""
fhist(; options...) = transformation(FHistAnalysis(; options...))
end
