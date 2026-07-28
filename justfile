default:
    just --list

install-julia-deps:
    #!/usr/bin/env -S julia --threads=auto --project=.
    using Pkg
    Pkg.instantiate()

figures:
    #!/usr/bin/env -S julia --threads=auto --project=.
    foreach(i -> include("scripts/fig$i.jl"), 1:5)

# Figure names must match the \includegraphics keys in overleaf/index.tex
cp-figures:
    cp figures/flux_with_fit.pdf overleaf/fig1-flux_with_fit.pdf
    cp figures/n_mlt_mlat.pdf overleaf/fig2-n_mlt_mlat.pdf
    cp figures/e30_flux_ratio_mlt_mlat_median.pdf overleaf/fig3-e30_flux_ratio_mlt_mlat_median.pdf
    cp figures/key_params_mlt_mlat_median.pdf overleaf/fig4-key_params_mlt_mlat_median.pdf
    cp figures/pdf_params_ae_mlt_mlat.pdf overleaf/fig5-pdf_params_ae_mlt_mlat.pdf

sparthb:
    quarto render sparthb.qmd --to pptx
    open _site/sparthb.pptx

flux-model-test:
    #!/usr/bin/env -S julia --threads=auto --project=lib/PrecipitatingFluxModels
    using Pkg
    Pkg.test()
