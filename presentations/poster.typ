#import "@preview/peace-of-posters:0.5.6" as pop
#import "@preview/unify:0.7.1": unit
#import "@preview/mitex:0.2.6": *

#set page("a0", margin: 1cm, flipped: true)
#pop.set-poster-layout(pop.layout-a0)
#pop.set-theme(pop.uni-fr)
#set text(size: pop.layout-a0.at("body-size"))
#let box-spacing = 1.2em
#set columns(gutter: box-spacing)
#set block(spacing: box-spacing)
#pop.update-poster-layout(spacing: box-spacing)

#pop.title-box(
  "Statistics of contribution of energetic electron precipitations: ELFIN and DMSP observations",
  // or "Peculiarities of precipitating electron spectra : DMPS & ELFIN combined dataset"
  authors: "Zijin Zhang, Anton V. Artemyev, Vassilis Angelopoulos",
  institutes: "University of California, Los Angeles",
  // keywords: "Peace, Dove, Poster, Science",
)



#columns(3, [
  #pop.column-box(heading: "Introduction - Energetic Electron Precipitation")[
    Energetic electron precipitation (EEP) couples the magnetosphere and ionosphere by creating ionization and controlling high‑latitude conductances. Empirical conductance models often rely on DMSP particle spectra as input. This work combines ELFIN precipitating electron measurements with DMSP spectra at close conjunction to extend the DMSP‑based precipitation picture toward higher energies.

    // Energetic electron precipitation couples magnetosphere and ionosphere
    // DMSP provides low-energy measurements (30 keV)
    // ELFIN extends to higher energies (50 keV to ~6 MeV)

    #figure(
      image("../figures/flux_with_fit.png", width: 100%),
      caption: [Example of a satellite conjunction event observed between 2021-12-01T22:19 and 2021-12-01T22:28. Panels (a–c) are plotted against magnetic latitude (MLAT): (a) the magnetic local time (MLT); (b) the precipitating electron flux measured by ELFIN; and (c) the electron flux measured by DMSP. Panels (d-e) display the precipitating flux spectra from both satellites, averaged over 0.5° MLAT bins.],
    )

    #mi(
      "j(E) = j_{EP} + j_κ = A_{EP} (\frac{E}{E_0})^{-γ}\,\exp(-\frac{E}{E_{EP}}) + A_κ \frac{E}{E_0} (1 + \frac{E}{κ E_κ})^{-κ-1}",
    )
  ]


  #pop.column-box()[
    #figure(
      image("../figures/n_mlt_mlat.png", width: 95%),
      caption: [
        The total number of MLAT-averaged spectral samples across different MLT and MLAT regions, sorted the AE index (2,754 conjunction events, with a total of 36,380 spectral fits).
      ],
    )
  ]

  // #pop.column-box(heading: "Model - Predict Precipitating Electron Flux")[

  //   Final products: 1. Event catalog and 2. empirical, data-driven model.

  //   #set text(30pt)

  //   ```julia
  //   using DEEEP # DmspElfinEnergeticElectronPrecipitation

  //   model = load_model()

  //   # Evaluate flux at specific geophysical conditions: Location: 65° MLat, 6 MLT (dawn sector), Activity: Moderate (AE = 150 nT)
  //   j_Efunc = model(; mlat=65.0, mlt=6.0, ae=150.0)

  //   E = 10.0 # [keV]
  //   flux_10keV = j_Efunc(E) # Calculate flux at a specific energy

  //   energies = 10 .^ range(log10(0.03), log10(1000), length=100)  # 0.03 - 1000 keV
  //   spectrum = j_Efunc.(energies) # Generate the full energy spectrum
  //   ```
  // ]

  #colbreak()

  #pop.column-box(heading: "Results - MLAT-MLT-AE Dependence")[

    #set text(36pt)

    Our analysis reveals that energetic flux, number flux, average energy and kappa parameter exhibit a strong dependence on MLAT, MLT, and geomagnetic activity.

    - #box(
        stroke: rgb("#8B4513"),
        fill: rgb("#FFE5B4"),
        inset: 8pt,
        radius: 4pt,
      )[The strongest energy and number fluxes occur in the *dawn* sector at high magnetic latitudes ($65°$–$75°$), the relative contribution of energetic particles exhibits an *anti-correlated* behavior.]

    - The energetic ($>30$ keV) component makes a substantial contribution to the total energy flux in the *dusk* sector, emphasizing the importance of accounting for this population to accurately characterize magnetosphere–ionosphere coupling.

    Building on these results, we developed an empirical, data-driven model to characterize precipitating electron fluxes from the thermal regime (∼ 1 keV, DMSP) to relativistic energies (∼ 1 MeV, ELFIN).
  ]

  #pop.column-box()[
    #figure(
      image("../figures/e30_flux_ratio_mlt_mlat_median.png", width: 100%),
      caption: [
        Median distributions of (a) total energy flux [$unit("keV cm^-2 s^-1 sr^-1")$], (b) total number flux [$unit("cm^-2 s^-1 sr^-1")$], and (c-d) fractional contribution of energetic particles (>30 keV) to total energy and number flux.
      ],
    )

    Increases in the AE index are accompanied by enhanced total energy and number fluxes, with the region of significant precipitation expanding toward lower magnetic latitudes and a broader MLT range.

    In contrast, the fractional contribution of the energetic component decreases during periods of high AE activity relative to quieter intervals.
  ]

  #colbreak()

  #pop.column-box()[

    #figure(
      image("../figures/key_params_mlt_mlat_median.pdf", width: 70%),
      caption: [Median distributions of (a) averaged energy ($overline(E)$) [$"keV"$] and (b) kappa parameter ($κ$) as functions of MLT, MLAT and AE index levels.],
    )

    #figure(
      image("../figures/pdf_params_ae_mlt_mlat.png", width: 100%),
      caption: [
        Probability density functions of (a) the energy flux ratio $J_E^{>30 "keV"}/J_E$, (b) the mean energy $overline(E)$ [keV], and (c) the kappa parameter $kappa$. Distributions are organized by three AE index ranges (0–100, 100–300, and $>300$ nT), two MLT sectors (12–21 and 0–9), and two MLAT bands ($60^∘$–$65^∘$ and $>65^∘$).
      ],
    )
  ]
])


// #pop.bottom-box()[
//   Peace of posters:
//   https://jonaspleyer.github.io/peace-of-posters/
// ]
