default:
    just --list

install-julia-deps:
    #!/usr/bin/env -S julia --threads=auto --project=.
    using Pkg
    Pkg.develop([
        PackageSpec("Speasy"),
        PackageSpec(url="https://github.com/Beforerr/MadrigalWeb.jl"),
        PackageSpec(url="https://github.com/Beforerr/GeoAACGM.jl"),
        PackageSpec(url="https://github.com/Beforerr/PySPEDAS.jl"),
    ])
    Pkg.instantiate()