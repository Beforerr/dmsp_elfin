## Source layout

- `src/fit.jl` — core fitting: `fit_two_step`, `init_guess`, `fit_two_flux`
- `src/sciml.jl` — NonlinearSolve (SciML) based solver (primary path)
- `src/model.jl`, `src/powerlaw.jl`, `src/kappa.jl` — model definitions
- `benchmark.jl` — benchmark with real event data (id=16, probe=a, 2022-07-02)

## NonlinearSolve performance notes

Benchmark target: `fit(TwoStepModel{PowerLawExpCutoff2, TransformKappaDistribution}, flux, energies)`
Current baseline (2026-05-07): **~4.0ms, 24.5k allocs, 1.6 MB** (no refine); **~7.3ms, 35k allocs** (refine).

Previous (before warm-start): ~7.2ms, 35k allocs (no refine); ~17.8ms, 59k allocs (refine).

### Warm-start across E_tran candidates

The outer `fit()` loop tries each energy channel as a candidate transition energy. Adjacent candidates share similar optimal M1/M2 solutions (one boundary point shifts). Passing the previous iteration's solution as `u0` cuts solver iterations.

### Choosing the right solver path

Three dispatch paths exist for `NonlinearLeastSquaresProblem`:

| Path | Speed | Allocs | Reason |
|---|---|---|---|
| IIP + `AutoSpecialize` | slow | low | `FunctionWrappersWrapper` forces chunk_size=1; 3 serial FD passes for 3 params |
| OOP + `AutoSpecialize` | fast | medium | No wrapper; chunk_size=3 (1 FD pass). But `NonlinearFunctionWrapper{false}` does `f(p,X) .- target` → 2 allocs per residual eval |
| **IIP + `FullSpecialize`** | **fast** | **low** | No wrapper → chunk_size=3. `DI.jacobian!` and `f!` both write in-place |

### Key implementation pattern

```julia
f! = (r, u, _E) -> (m = _eval_model(Model(u)); @. r = log_eval(m, _E) - log_y; nothing)
nf = NonlinearFunction{true, SciMLBase.FullSpecialize}(f!; resid_prototype = similar(log_y))
prob = NonlinearLeastSquaresProblem(nf, u0, E)
```

- `resid_prototype` is **required** when `length(residual) ≠ length(u0)` (overdetermined NLLS); without it, solver allocates `fu` as `similar(u0)` (wrong size).