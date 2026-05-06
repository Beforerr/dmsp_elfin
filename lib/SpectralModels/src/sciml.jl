using NonlinearSolveFirstOrder: FastShortcutNLLSPolyalg

# https://docs.sciml.ai/NonlinearSolve/stable/solvers/nonlinear_least_squares_solvers/
#
function sciml_log_fit(Model, E, y; alg = FastShortcutNLLSPolyalg(), u0 = init_guess(Model, E, y), kw...)
    log_y = log.(y)
    f! = (r, u, _E) -> (m = _eval_model(Model(u)); @. r = log_eval(m, _E) - log_y)
    nf = NonlinearFunction{true, _FS}(f!; resid_prototype = similar(log_y))
    prob = NonlinearLeastSquaresProblem(nf, u0, E)
    sol = solve(prob, alg; kw...)
    return Model(sol.u)
end
