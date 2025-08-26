# Example demo

```{julia}
using DmspElfinConjunction
using DmspElfinConjunction
using DMSP, ELFIN
using SPEDAS
using SpacePhysicsMakie
using Dates, TimeseriesUtilities
using DataFrames, DataFramesMeta, DimensionalData

using SpacePhysicsMakie: set_if_valid!
import DmspElfinConjunction.YLabel as 𝒀

# Use CairoMakie in CI, GLMakie otherwise
if haskey(ENV, "CI") && ENV["CI"] == "true"
    using CairoMakie
    CairoMakie.activate!()
else
    using GLMakie
    GLMakie.activate!()
end
```

```@example demo
probe = "a"
trange = ["2022-07-02", "2022-07-03"]

# Load ELFIN EPD data
elx_flux = ELFIN.precipitating_flux(trange, probe)
elx_gei = ELFIN.gei(trange, probe)
elx_gei = tinterp(elx_gei, elx_flux.anti)
elx_aacgm = gei2aacgm(elx_gei)
elx_mlat = elx_aacgm.mlat

dmsp_flux = DMSP.load(trange, 16, "el_d_flux")
dmsp_geod = DMSP.geod(trange, 16)
dmsp_aacgm = geod2aacgm(dmsp_geod)
dmsp_mlat = dmsp_aacgm.mlat

# use keV as the basic unit for energy dimension
dmsp_flux = set(dmsp_flux, Y => dims(dmsp_flux, Y).val .* 1e-3)
# set the flux unit to 1/cm²/s/sr/MeV
dmsp_flux *= 1e6
set_if_valid!(dmsp_flux.metadata, :yunit => "keV", :description => "Diff electron num flux (1/cm^2/s/sr/MeV)")
replace!(dmsp_flux, 0 => NaN)
```

```@example demo
for x in (elx_flux.para, elx_flux.anti, elx_flux.prec, dmsp_flux)
    set_if_valid!(x.metadata, :yscale => log10, :ylabel => "Energy (keV)", :colorscale => log10, :colorrange => (1e3, 1e11))
end

elx_mlat = SPEDAS.setmeta(elx_mlat, :label => "ELFIN")
set_if_valid!(dmsp_mlat.metadata, :labels => "DMSP")

trial_tr = ["2022-07-02T04:50", "2022-07-02T05:20"]
zoom_tr = timerange(tview(elx_flux.para, trial_tr))
tvars = (
    elx_flux.para, elx_flux.anti, elx_flux.prec,
    dmsp_flux,
    [elx_mlat, dmsp_mlat],
)
tplot(tvars, zoom_tr...)
```

```@example demo
elx_df = get_flux_by_mlat(elx_flux.prec, elx_mlat, zoom_tr)
dmsp_df = get_flux_by_mlat(dmsp_flux, dmsp_mlat, trial_tr)

join_df = leftjoin(dmsp_df, elx_df, on=:mlat, makeunique=true)
sort!(join_df, :mlat)
dropmissing!(join_df)
```

Fit the two-step model to the flux data on a single MLAT value.

```@example demo
mlats = [-61, -63, -66, -69]
flux_threshold = 100

sdf = @rsubset(join_df, :mlat ∈ mlats)
idx = 3
flux1 = sdf[idx, :].flux
flux2 = sdf[idx, :].flux_1
modelType = TwoStepModel{PowerLawExpCutoff, KappaDistribution}
res = fit_two_flux(modelType, flux1, flux2; flux_threshold)

f = Figure()
ax = plot_spectra(f[1, 1], flux1, flux2, res.model)
ylims!(ax, 1.0e1, 1.0e12)
axislegend(ax; position=:lb)
hlines!(ax, flux_threshold; color = :black, linestyle = :dash)
ax.title = "Two-Step Model Fit (MLAT = $(round(sdf[1, :].mlat, digits=1))°)"
f
```

## Fit all MLT values and analyze parameter variation

```@example demo
# Function to fit flux for a single row and extract parameters
# Fit all rows in the dataframe
@info "Fitting $(nrow(join_df)) MLT values..."

flux_threshold = 200
@rtransform! join_df $AsTable = fit_two_flux(:flux, :flux_1; flux_threshold)
successful_fits = filter(r -> r.success, join_df; view=true)
```


```@example demo
# Make MLAT as the x axis
f = Figure()
colorrange = (1e3, 1e11)
p1 = plot_flux_by_mlat(f[1, 1], elx_flux.prec, elx_mlat, zoom_tr; colorrange)
p2 = plot_flux_by_mlat(f[2, 1], dmsp_flux, dmsp_mlat, trial_tr; colorrange)
xlims!(p1.axis, -54, -75)
xlims!(p2.axis, -54, -75)
Colorbar(f[1:2, 2], p1.plot; label=𝒀.nflux)

for ax in (p1.axis, p2.axis)
    vlines!(ax, mlats)
end

sdf = @rsubset(join_df, :mlat ∈ mlats)
axs = plot_spectra(f[3, 1:end], sdf)

hlines!.(axs, flux_threshold; color = :black, linestyle = :dash)
axislegend.(axs; position=:lb)
hidexdecorations!(p1.axis; grid=false)
f
```

## Comprehensive Parameter Analysis

Plotting model-calculated flux vs MLAT and how the PowerLawExpCutoff parameters change over MLAT.

1. Spatial Variations: How particle precipitation characteristics change with magnetic latitude
2. Spectral Hardness Trends: Whether spectra become harder/softer at different MLATs
3. Cutoff Energy Patterns: How the high-energy cutoff varies spatially
<!-- 4. Parameter Correlations: Relationships between spectral parameters that reveal physical processes -->

```@example demo
modeled_fluxes = @with successful_fits begin
    flux_modeled = stack(:model) do model
        model.(energies)
    end
    DimArray(flux_modeled', (X(:mlat), Y(energies)))
end
```


```@example demo
# Make MLAT as the x axis
f = Figure(; size=(1200, 1000))
p1 = plot_flux_by_mlat(f[1, 1:2], elx_flux.prec, elx_mlat, zoom_tr; colorrange)
p2 = plot_flux_by_mlat(f[2, 1:2], dmsp_flux, dmsp_mlat, trial_tr; colorrange)
p3 = plot_flux_by_mlat(f[3:4, 1:2], modeled_fluxes; colorrange)
Colorbar(f[1:4, 3], p3.plot; label=𝒀.nflux)

for ax in (p1.axis, p2.axis, p3.axis)
    xlims!(ax, -54, -75)
end
hidexdecorations!.((p1.axis, p2.axis); grid=false)

@with successful_fits begin
    plot_parameters_variation(f[1:4, 4:5], :mlat, :model, :n_points; scores = :score)
end

f
```