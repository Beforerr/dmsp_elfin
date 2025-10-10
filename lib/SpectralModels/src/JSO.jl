using ADNLPModels, CaNNOLeS, JSOSolvers

# JSO (JuliaSmoothOptimizers) interface for nonlinear least squares
"""
    jso_nls_fit(Model, E, y; solver=cannoles, bounds=nothing, kw...)

Fit spectral model using JSO nonlinear least squares solvers.

Uses ADNLPModels.jl for problem modeling and CaNNOLeS.jl for solving.

# Example
```julia
using CaNNOLeS
model = jso_nls_fit(KappaDistribution, energies, flux)

# Or with JSOSolvers
using JSOSolvers
model = jso_nls_fit(KappaDistribution, energies, flux, solver=trunk)
```

# References
- https://jso.dev/ADNLPModels.jl/dev/
- https://jso.dev/JSOSolvers.jl/latest/solvers/
"""
function jso_nls_fit(Model, E, y; bounds = nothing, kw...)
    # Residual function for nonlinear least squares
    observed = nm.log.(y)
    function F(x)
        m = Model(x)
        return observed .- log_eval.(m, E)
    end

    p0 = init_guess(Model, E, y)

    # Create constrained ADNLPModel
    nls = if bounds !== nothing
        lvar, uvar = bounds
        ADNLSModel(F, p0, length(E), lvar, uvar)
    else
        ADNLSModel(F, p0, length(E))
    end

    # Solve using JSO solver
    stats = if isnothing(bounds)
        solver = TrunkSolverNLS(nls)
        JSOSolvers.solve!(solver, nls; kw...)
    else
        tron(nls)
    end

    return Model(stats.solution)
end