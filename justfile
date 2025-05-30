default:
    just --list

install-julia-deps:
    #!/usr/bin/env -S julia --threads=auto --project=.
    using Pkg
    Pkg.develop([
        PackageSpec("Speasy.jl"),
        PackageSpec(path="./MadrigalWeb.jl"),
        PackageSpec(url="https://github.com/Beforerr/PySPEDAS.jl"),
    ])
    Pkg.instantiate()