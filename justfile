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

sparthb:
    quarto render sparthb.qmd --to pptx
    quarto render article.qmd --to pptx
    open _site/sparthb.pptx

cp-figures:
    cp figures/flux_with_fit.pdf overleaf/flux_with_fit.pdf
    cp figures/n_mlt_mlat.pdf overleaf/n_mlt_mlat.pdf
    cp figures/e30_flux_ratio_mlt_mlat_median.pdf overleaf/e30_flux_ratio_mlt_mlat_median.pdf
    cp figures/key_params_mlt_mlat_median.pdf overleaf/key_params_mlt_mlat_median.pdf
    cp figures/pdf_params_ae_mlt_mlat.pdf overleaf/pdf_params_ae_mlt_mlat.pdf

# PrecipitatingFluxModels package tasks
flux-model-create:
    #!/usr/bin/env -S julia --threads=auto --project=docs
    include("scripts/create_flux_model.jl")

flux-model-test:
    #!/usr/bin/env -S julia --threads=auto --project=lib/PrecipitatingFluxModels
    using Pkg
    Pkg.test()

flux-model-examples:
    #!/usr/bin/env -S julia --threads=auto --project=lib/PrecipitatingFluxModels
    include("lib/PrecipitatingFluxModels/examples/basic_usage.jl")

