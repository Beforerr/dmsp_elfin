mlt2x(mlt) = (mlt + 12) % 24
mlt2rad(x) = x / 24 * 2π
mlt2rad(x::AbstractString) = parse(Float64, String(x)) |> mlt2rad
rad2mlt(x) = x / 2π * 24

const mlt_tickformat = x -> @. string(round(Int, (x + 12) % 24))

function polarplot!(ax, mlts, mlats, values; kw...)
    θ = mlt2rad.(mlts)
    r = @. 90 - _mlat(mlats)
    # Create matrix with NaN for missing combinations
    unique_θ = sort(unique(θ))
    unique_r = sort(unique(r))
    matrix = fill(NaN, length(unique_θ), length(unique_r))
    for (t, rad, val) in zip(θ, r, values)
        i = findfirst(==(t), unique_θ)
        j = findfirst(==(rad), unique_r)
        if !isnothing(i) && !isnothing(j)
            matrix[i, j] = val
        end
    end

    return voronoiplot!(ax, unique_θ, unique_r, matrix; show_generators = false, strokewidth = 0, kw...)
end

function polarplot(fp, args...; kw...)
    ax = PolarAxis(
        fp;
        rtickformat = x -> string.(90 .- x),
        thetaticks = AngularTicks(24 / 2π, "")
    )
    return ax, polarplot!(ax, args...; kw...)
end

# Reorder MLT levels: 12→24, then 0→12
function mltmlatplot!(ax, mlts, mlats, values; kw...)
    x = mlt2x.(mlts)
    y = _mlat.(mlats)
    return heatmap!(ax, x, y, values; kw...)
end

function mltmlatplot(fp, args...; colorscale = identity, axis = (;), kw...)
    ax = Axis(fp; xlabel = "MLT", ylabel = "MLAT", xtickformat = mlt_tickformat, axis...)
    return ax, mltmlatplot!(ax, args...; colorscale, kw...)
end
