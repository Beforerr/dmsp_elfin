# Example demo

```{julia}
using DmspElfinConjunction
using DmspElfinConjunction.ELFIN
using DMSP
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
@info dmsp_flux.dims[2].val.metadata
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
    set_if_valid!(x.metadata,:yscale => log10, :ylabel => "Energy (keV)", :colorscale => log10, :colorrange => (1e3, 1e11))
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

```@example demo
# Make MLAT as the x axis
f = Figure()
colorrange = (1e3, 1e11)
p1 = plot_flux_by_mlat(f[1, 1], elx_flux.prec, elx_mlat, zoom_tr; colorrange)
p2 = plot_flux_by_mlat(f[2, 1], dmsp_flux, dmsp_mlat, trial_tr; colorrange)
xlims!(p1.axis, -54, -75)
xlims!(p2.axis, -54, -75)
Colorbar(f[1:2, 2], p1.plot; label=𝒀.nflux)

mlats = [-61, -63, -66, -69]

for ax in (p1.axis, p2.axis)
    vlines!(ax, mlats)
end

sdf = @rsubset(join_df, :mlat ∈ mlats)
plot_spectra(f[3, 1:end], sdf)

hidexdecorations!(p1.axis; grid=false)
f
```

```@example demo
flux1 = sdf[1, :].flux
flux2 = sdf[1, :].flux_1
ff = vcat(flux1.data, flux2.data)
ee = vcat(flux1.dims[1].val, flux2.dims[1].val)
energies = ee
ff, ee = remove_nan(ff, ee)

model, flux_modeled, Emin, score = fit_flux_two_step(ff, ee)

f = Figure()
ax = plot_spectra(f[1, 1], flux1, flux2, model)
ylims!(ax, 1.0e1, 1.0e11)
axislegend(ax; position=:lb)
ax.title = "Two-Step Model Fit (MLAT = $(round(sdf[1, :].mlat, digits=1))°)"
f
```

## Fit all MLT values and analyze parameter variation

```@example demo
# Function to fit flux for a single row and extract parameters
# Fit all rows in the dataframe
@info "Fitting $(nrow(join_df)) MLT values..."

@rtransform! join_df $AsTable = fit_row_parameters(:flux, :flux_1)
successful_fits = filter(r -> r.success, join_df; view=true)
```

## Comprehensive Parameter Analysis

Plotting model-calculated flux vs MLAT and how the PowerLawExp parameters change over MLAT.

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