# DMSP-ELFIN Conjunction Study

[![Paper](https://img.shields.io/badge/JGR%20Space%20Physics-10.1029%2F2026JA035299-blue)](https://doi.org/10.1029/2026JA035299)
[![Software](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.20152212-1682D4)](https://doi.org/10.5281/zenodo.20152212)
[![Dataset](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.20132654-1682D4)](https://doi.org/10.5281/zenodo.20132654)
[![Julia](https://img.shields.io/badge/-Julia-9558B2?logo=julia&logoColor=white)](https://julialang.org/)

Analysis code and manuscript sources for:

> Zhang, Z., Artemyev, A. V., Angelopoulos, V., & Zhang, X.-J. (2026). Peculiarities of precipitating electron spectra revealed through a combined DMSP-ELFIN data set. *Journal of Geophysical Research: Space Physics*, 131(6), e2026JA035299. <https://doi.org/10.1029/2026JA035299>

Existing precipitation datasets rest on DMSP measurements below 30 keV and so miss the energetic electrons that dominate lower-ionosphere dynamics. This study pairs conjugate DMSP and ELFIN measurements to build a combined dataset of ~3.6 × 10⁴ electron spectra, fits each spectrum with an analytical composite model, and packages the result as an empirical model of precipitating flux as a function of MLAT, MLT, and AE.

## Products

| Product | Location |
| --- | --- |
| Event catalog (conjunctions, fluxes, fitted spectral parameters) | [Zenodo 10.5281/zenodo.20132654](https://doi.org/10.5281/zenodo.20132654) |
| DEEEP flux model (Julia + Python) | [`lib/PrecipitatingFluxModels`](./lib/PrecipitatingFluxModels/README.md) |
