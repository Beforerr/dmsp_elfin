default:
    just --list

install-julia-deps:
    #!/usr/bin/env -S julia --threads=auto --project=.
    using Pkg
    # Pkg.update()
    Pkg.develop([
        PackageSpec("Madrigal"),
        PackageSpec("TimeseriesUtilities"),
        PackageSpec(url="https://github.com/jishnub/SphericalHarmonics.jl"),
        PackageSpec(url="https://github.com/SciML/CurveFit.jl"),
    ])
    Pkg.instantiate()