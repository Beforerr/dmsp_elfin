using LsqFit

# https://github.com/JuliaNLSolvers/LsqFit.jl
function LsqFit_log_fit(Model, E, y; kw...)
    f = (E, p) -> (m = Model(p); log_eval.(m, E))
    p0 = init_guess(Model, E, y)
    log_y = log.(y)
    try
        fit = LsqFit.curve_fit(f, E, log_y, p0; kw...)
        return Model(fit.param)
    catch e
        @info log_y p0
        throw(e)
    end
end
