#import "@preview/touying:0.6.1": *
#import "@preview/unify:0.7.1": unit
#import themes.university: *

#show: university-theme.with(
  aspect-ratio: "16-9",
  config-info(
    title: ["Peculiarities of precipitating electron spectra: DMSP & ELFIN combined dataset"],
  ),
)
// #show: simple-theme.with(aspect-ratio: "16-9")

== Introduction

#slide(composer: (1fr, 2fr))[
  #set text(20pt)

  - Energetic electron precipitation couples magnetosphere and ionosphere
  - DMSP provides low-energy measurements (30 keV)
  - ELFIN extends to higher energies (50 keV to ~6 MeV)
][
  #show figure.caption: set text(15pt)
  #move(
    dy: -1.5cm,
    figure(
      image("../figures/flux_with_fit.png"),
      caption: [Example of a satellite conjunction event observed between 2021-12-01T22:19 and 2021-12-01T22:28. Panels (a–c) are plotted against magnetic latitude (MLAT): (a) the magnetic local time (MLT); (b) the precipitating electron flux measured by ELFIN; and (c) the electron flux measured by DMSP. Panels (d-e) display the precipitating flux spectra from both satellites, averaged over 0.5° MLAT bins.],
    ),
  )
]

== Results

#figure(
  image("../figures/n_mlt_mlat.png"),
  caption: [The total number of MLAT-averaged spectral data points across different MLT and MLAT regions, separated by varying levels of the AE index.],
)

2,754 conjunction events identified.

#slide(composer: (4fr, 3fr))[

  #show figure.caption: set text(13pt)

  #figure(
    image("../figures/e30_flux_ratio_mlt_mlat_median.png"),
    caption: [Median distributions of (a) total energy flux [$unit("keV cm^-2 s^-1 sr^-1")$], (b) total number flux [$unit("cm^-2 s^-1 sr^-1")$], and (c-d) fractional contribution of energetic particles (>30 keV) to total energy and number flux.],
  )
][
  #show figure.caption: set text(13pt)
  #pause
  #figure(
    image("../figures/key_params_mlt_mlat_median.pdf", width: 90%),
    caption: [Median distributions of (a) averaged energy ($overline(E)$) [$"keV"$] and (b) kappa parameter ($κ$) as functions of MLT and MLAT, sorted by AE index levels.],
  )
]

#show figure.caption: set text(15pt)

#figure(
  image("../figures/pdf_params_ae_mlt_mlat.png"),
  caption: [Probability density functions for (a) energy flux ratio $J_E^(>30 "keV")/J_E$, (b) mean energy $overline(E)$ and (c) kappa parameter ($κ$). Data organized by AE ranges, MLT sectors, and MLAT bands.],
)
