function remove_nan(Xs...)
    idx = mapreduce(.&, Xs) do x
        .!isnan.(x)
    end
    return map(x -> x[idx], Xs)
end

energies(x) = parent(x.dims[1].val)


function dist(f, x, y)
    d = Inf
    for xi in x
        for yi in y
            d = min(d, f(xi, yi))
        end
    end
    return d
end

function dist(x, y)
    return dist(abs ∘ (-), x, y)
end

function mlt_dist(x, y)
    d = abs(x - y)
    return min(d, 24 - d)
end
