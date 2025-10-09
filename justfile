default:
    just --list

install-julia-deps:
    #!/usr/bin/env -S julia --threads=auto --project=.
    using Pkg
    # Pkg.update()
    Pkg.develop([
        PackageSpec("Madrigal"),
        PackageSpec("TimeseriesUtilities"),
    ])
    Pkg.instantiate()