using NonlinearSolve
using CurveFit

function loss_function!(resid, u, p)
    Model, E, log_y = p
    m = Model(u)
    resid .= log_y .- log_eval.(m, E)
    return
end

function loss_function(u, p)
    m = p.Model(u)
    return p.log_y .- log_eval.(m, p.E)
end

# https://docs.sciml.ai/NonlinearSolve/stable/solvers/nonlinear_least_squares_solvers/
function sciml_log_fit(Model, E, y; alg = FastShortcutNLLSPolyalg(), kw...)
    u0 = init_guess(Model, E, y)
    log_y = log.(y)
    _prob = :NonlinearCurveFitProblem
    if _prob == :NonlinearLeastSquaresProblem # slower but allocations are smaller
        p = (Model, E, log_y)
        prob = NonlinearLeastSquaresProblem(NonlinearFunction(loss_function!, resid_prototype = similar(y)), u0, p)
    elseif _prob == :NonlinearCurveFitProblem
        f = (u, E) -> (m = Model(u); log_eval.(m, E))
        prob = NonlinearCurveFitProblem(f, u0, E, log_y)
    end
    sol = solve(prob, alg; kw...)
    return Model(sol.u)
end
