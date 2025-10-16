### PrecipitatingFluxModels.jl Usage Guide

Complete guide for creating and using empirical flux models.

## Table of Contents

1. [For Model Creators](#for-model-creators)
2. [For End Users](#for-end-users)
3. [Advanced Usage](#advanced-usage)
4. [API Reference](#api-reference)

## For End Users

### Computing Energy Spectra

```julia
# Plot
using CairoMakie
fig, ax, ln = lines(energies, spectrum;
    axis=(xscale=log10, yscale=log10,
          xlabel="Energy (keV)",
          ylabel="Flux (cm⁻² s⁻¹ sr⁻¹ keV⁻¹)")
)
```

### Uncertainty Quantification

```julia
# Get median and percentiles
flux_med = evaluate_flux(model, 10.0; mlat=65.0, mlt=6.0, ae=150.0, stat=:median)
flux_low = evaluate_flux(model, 10.0; mlat=65.0, mlt=6.0, ae=150.0, stat=:q25)
flux_high = evaluate_flux(model, 10.0; mlat=65.0, mlt=6.0, ae=150.0, stat=:q75)

println("Flux: $(flux_med) [$(flux_low) - $(flux_high)]")
```

### Integrated Fluxes

```julia
# Number flux (particles cm⁻² s⁻¹ sr⁻¹)
J = number_flux(model, 0.03, 100.0;  # 0.03-100 keV
    mlat=65.0, mlt=6.0, ae=150.0
)

# Energy flux (keV cm⁻² s⁻¹ sr⁻¹)
JE = energy_flux(model, 0.03, 100.0;
    mlat=65.0, mlt=6.0, ae=150.0
)

println("Number flux: $J")
println("Energy flux: $JE")
```

---

## Advanced Usage

### Spatial Mapping

Create maps of flux across MLat and MLT:

```julia
using CairoMakie

mlats = 60:0.5:75
mlts = 0:0.5:23.5
energy = 10.0  # keV
ae = 150.0     # nT

# Compute flux map
flux_map = [evaluate_flux(model, energy; mlat=mlat, mlt=mlt, ae=ae)
            for mlat in mlats, mlt in mlts]

# Plot
fig, ax, hm = heatmap(mlts, mlats, flux_map';
    colorscale=log10,
    colormap=:viridis,
    axis=(xlabel="MLT (hours)", ylabel="MLat (degrees)",
          title="$energy keV Flux (AE=$ae nT)")
)
Colorbar(fig[1,2], hm, label="Flux (cm⁻² s⁻¹ sr⁻¹ keV⁻¹)")
save("flux_map.png", fig)
```

### Activity Dependence

Compare flux levels across different geomagnetic conditions:

```julia
ae_levels = [50.0, 150.0, 400.0]  # Quiet, moderate, active
energy = 10.0
mlat, mlt = 65.0, 6.0

fig = Figure()
ax = Axis(fig[1,1];
    xlabel="AE (nT)",
    ylabel="Flux (cm⁻² s⁻¹ sr⁻¹ keV⁻¹)",
    yscale=log10
)

fluxes = [evaluate_flux(model, energy; mlat=mlat, mlt=mlt, ae=ae)
          for ae in ae_levels]

scatter!(ax, ae_levels, fluxes; markersize=15)
lines!(ax, ae_levels, fluxes)

save("ae_dependence.png", fig)
```

### Spectral Components

Examine individual spectral contributions:

```julia
components = flux_components(model, 50.0;
    mlat=65.0, mlt=6.0, ae=150.0
)

println("Low-energy (ExpPow): $(components.low_energy)")
println("High-energy (Kappa): $(components.high_energy)")
println("Total: $(components.total)")
println("Kappa fraction: $(components.high_energy/components.total * 100)%")
```

### Accessing Raw Parameters

For custom analysis, access model parameters directly:

```julia
params = flux_parameters(model; mlat=65.0, mlt=6.0, ae=150.0)

println("κ: $(params.κ)")
println("E_c1: $(params.E_c1) keV")
println("E_c2: $(params.E_c2) keV")
println("γ: $(params.γ)")
println("Transition energy: $(params.Emin) keV")
println("Samples: $(params.n_samples)")

# Convert to SpectralModel for custom evaluation
spectral_model = to_spectral_model(params)
custom_flux = spectral_model(25.0)  # Evaluate at 25 keV
```

### Complete Example with Plot

```julia
using PrecipitatingFluxModels
using CairoMakie

# Create energy grid
energies = 10 .^ range(log10(0.03), log10(1000), length=100)

# Evaluate spectrum with uncertainty
spectrum_med = flux_spectrum(model, energies; mlat=65.0, mlt=6.0, ae=150.0, stat=:median)
spectrum_low = flux_spectrum(model, energies; mlat=65.0, mlt=6.0, ae=150.0, stat=:q25)
spectrum_hi = flux_spectrum(model, energies; mlat=65.0, mlt=6.0, ae=150.0, stat=:q75)

# Plot
fig = Figure(size=(800, 600))
ax = Axis(fig[1,1];
    xlabel="Energy (keV)",
    ylabel="Differential Flux (cm⁻² s⁻¹ sr⁻¹ keV⁻¹)",
    xscale=log10, yscale=log10,
    title="Precipitating Electron Flux (MLat=65°, MLT=6h, AE=150 nT)"
)

# Plot median with uncertainty band
band!(ax, energies, spectrum_low, spectrum_hi; color=(:steelblue, 0.3), label="25-75%")
lines!(ax, energies, spectrum_med; color=:steelblue, linewidth=2, label="Median")

axislegend(ax; position=:lb)
save("flux_spectrum.png", fig)
```

---

## API Reference

### Flux Evaluation

```julia
evaluate_flux(model, energy; mlat, mlt, ae, stat=:median)
    # Single energy evaluation
    # Returns: Float64 (cm⁻² s⁻¹ sr⁻¹ keV⁻¹)

flux_spectrum(model, energies; mlat, mlt, ae, stat=:median)
    # Energy spectrum
    # Returns: Vector{Float64} (cm⁻² s⁻¹ sr⁻¹ keV⁻¹)

number_flux(model, E_min, E_max; mlat, mlt, ae, stat=:median)
    # Integrated number flux
    # Returns: Float64 (cm⁻² s⁻¹ sr⁻¹)

energy_flux(model, E_min, E_max; mlat, mlt, ae, stat=:median)
    # Integrated energy flux
    # Returns: Float64 (keV cm⁻² s⁻¹ sr⁻¹)

flux_components(model, energy; mlat, mlt, ae, stat=:median)
    # Individual components
    # Returns: NamedTuple(low_energy, high_energy, total)
```

### Parameters

```julia
flux_parameters(model; mlat, mlt, ae, stat=:median)
    # Get FluxParameters object
    # Returns: FluxParameters
```

### Common Arguments

- `mlat`: Magnetic latitude (degrees, -90 to 90)
- `mlt`: Magnetic local time (hours, 0-24)
- `ae`: AE index (nT, typically 0-1000+)
- `stat`: Statistics to use (`:median`, `:q25`, `:q75`)

---

## Tips and Best Practices

### For Model Creators

1. **Sufficient Statistics**: Ensure each spatial bin has enough observations (>20 recommended)
2. **Quality Control**: Filter out poor fits before computing statistics

### For End Users

1. **Interpolation**: Model uses nearest-neighbor interpolation; results most reliable near bin centers
2. **Energy Range**: Respect the model's valid energy range (check `model.energy_range`)
3. **Uncertainty**: Always consider using percentiles for uncertainty estimation
4. **Validation**: Compare model predictions with observations when available

**Out of bounds values:**
- Model will use nearest neighbor for out-of-range coordinates
- Extrapolation may not be reliable

**Missing columns error:**
- Ensure all required parameters have _median, _q25, _q75 columns
- Check column names match exactly