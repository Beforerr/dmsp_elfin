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

# the mean of [1,23] should be 0
cos24(x) = cos(x * 2π / 24)
sin24(x) = sin(x * 2π / 24)

function local_mlt_mean(x)
    real = mean(cos24, x)
    imag = mean(sin24, x)
    return mod(angle(real + imag * im) * 24 / 2π, 24)
end
