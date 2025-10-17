# PrecipitatingFluxModels.jl

A Julia package for modeling precipitating electron flux using empirical statistical models derived from DMSP and ELFIN satellite observations.

## Features

- 🌍 **Spatial Coverage**: Models flux as a function of magnetic latitude (MLat), magnetic local time (MLT), and geomagnetic activity (AE index)
- 📊 **Statistical Uncertainty**: Provides median values plus 25th and 75th percentiles for all parameters
- 🔬 **Two-Component Spectra**: Combines exponential power-law (thermal) and kappa distribution (suprathermal) models
- ⚡ **Fast Evaluation**: Pre-computed statistical parameters for quick flux calculations
- 💾 **Serializable**: Save and load models in JLD2 format

## Installation

```julia
using Pkg
Pkg.develop(path="path/to/PrecipitatingFluxModels")
```

## Quick Start

```julia
using DEEEP # DmspElfinEnergeticElectronPrecipitation

# Load the default model
model = load_flux_model()

# Evaluate flux at specific conditions
# Location: 65° MLat, 6 MLT (dawn sector), moderate activity (AE=150 nT)
flux_10keV = evaluate_flux(model, 10.0; mlat=65.0, mlt=6.0, ae=150.0)

# Get full energy spectrum
energies = 10 .^ range(log10(0.03), log10(1000), length=100)  # 0.03-1000 keV
spectrum = flux_spectrum(model, energies; mlat=65.0, mlt=6.0, ae=150.0)

# Compute integrated fluxes
Emin = 0.03
Emax = 100.0
J = number_flux(model, Emin, Emax; mlat=65.0, mlt=6.0, ae=150.0)  # cm⁻² s⁻¹ sr⁻¹
JE = energy_flux(model, Emin, Emax; mlat=65.0, mlt=6.0, ae=150.0)  # keV cm⁻² s⁻¹ sr⁻¹
```

## Usage Examples

### Basic Flux Evaluation

```julia
using PrecipitatingFluxModels

model = load_flux_model()

# Auroral zone, morning sector, quiet conditions
flux = evaluate_flux(model, 50.0; mlat=67.0, mlt=8.0, ae=50.0)
println("Flux at 50 keV: $flux cm⁻² s⁻¹ sr⁻¹ keV⁻¹")
```

### Uncertainty Quantification

```julia
# Get median and uncertainty bounds
flux_median = evaluate_flux(model, 50.0; mlat=67.0, mlt=8.0, ae=50.0, stat=:median)
flux_lower = evaluate_flux(model, 50.0; mlat=67.0, mlt=8.0, ae=50.0, stat=:q25)
flux_upper = evaluate_flux(model, 50.0; mlat=67.0, mlt=8.0, ae=50.0, stat=:q75)

println("Flux: $flux_median [$flux_lower - $flux_upper]")
```

### Spectral Components

```julia
# Examine contributions from different populations
components = flux_components(model, 50.0; mlat=67.0, mlt=8.0, ae=150.0)

println("Low-energy (thermal): ", components.low_energy)
println("High-energy (tail): ", components.high_energy)
println("Total: ", components.total)
```

### Spatial Mapping

```julia
using CairoMakie

# Create a spatial map of flux at 10 keV for moderate activity
mlats = 60:1:75
mlts = 0:1:23

flux_map = [evaluate_flux(model, 10.0; mlat=mlat, mlt=mlt, ae=150.0)
            for mlat in mlats, mlt in mlts]

fig, ax, hm = heatmap(mlts, mlats, flux_map';
    axis=(xlabel="MLT (hours)", ylabel="MLat (degrees)",
          title="10 keV Flux (AE=150 nT)"),
    colorscale=log10
)
Colorbar(fig[1,2], hm, label="Flux (cm⁻² s⁻¹ sr⁻¹ keV⁻¹)")
fig
```

## Model Parameters

The model uses a two-component spectral representation:

### Low Energy (0.03-~30 keV): Exponential Power Law
```
f₁(E) = A₁ · E^(-γ) · exp(-(E/E_c1))
```

### High Energy (~30-1000 keV): Kappa Distribution
```
f₂(E) = A₂ · [1 + E/(κ·E_c2)]^(-κ-1)
```

### Combined Model
```
f(E) = f₁(E) + f₂(E)
```

## Creating Your Own Model

If you have statistical analysis results, you can create a custom model:

```julia
using DataFrames, DataFramesMeta
using PrecipitatingFluxModels

# Prepare your statistics (see format requirements in docs)
stats_df = compute_your_statistics()  # Your analysis pipeline

# Create model
model = create_model_from_stats(stats_df;
    mlat_bins=50:1:80,
    mlt_bins=0:2:22,
    ae_bins=["[0, 100)", "[100, 300)", "[300, Inf)"],
    energy_range=(0.03, 1000.0),
    description="My Custom Model",
    version="1.0.0"
)

# Save for distribution
save_flux_model("my_model.jld2", model)

# Load later
my_model = load_flux_model("my_model.jld2")
```

## Data Format

Statistical parameter DataFrame must include for each spatial bin (MLat, MLT, AE):

**Spatial Binning:**
- `mlat_bin`: Magnetic latitude bin center (degrees)
- `mlt_bin`: Magnetic local time bin center (hours)
- `ae_bin`: AE index bin (categorical string)
- `n_samples`: Number of observations

**Model Parameters** (with `_median`, `_q25`, `_q75` suffixes):
- `κ`: Kappa parameter
- `log_A1`, `E_c1`, `γ`: ExpPow model parameters
- `log_A2`, `E_c2`: Kappa model parameters
- `Emin`: Transition energy (keV)
- `J1`, `J2`: Number flux components (cm⁻² s⁻¹ sr⁻¹)
- `JE1`, `JE2`: Energy flux components (keV cm⁻² s⁻¹ sr⁻¹)

## API Reference

### Core Functions

- `load_flux_model([filename])`: Load model from file or use default
- `evaluate_flux(model, energy; mlat, mlt, ae, stat)`: Evaluate at single energy
- `flux_spectrum(model, energies; mlat, mlt, ae, stat)`: Compute spectrum
- `number_flux(model, E_min, E_max; kwargs...)`: Integrated number flux
- `energy_flux(model, E_min, E_max; kwargs...)`: Integrated energy flux
- `flux_parameters(model; mlat, mlt, ae, stat)`: Get raw parameters
- `flux_components(model, energy; kwargs...)`: Get spectral components

### Model Creation

- `create_model_from_stats(df; kwargs...)`: Create from statistics
- `save_flux_model(filename, model)`: Save to file
- `get_ae_bin(ae_value, ae_bins)`: Find appropriate AE bin

## Citation

If you use this package in your research, please cite:

```
Zhang, Z. et al. (2024). Empirical Model of Precipitating Electron Flux
from DMSP and ELFIN Observations. [Journal/Preprint info]
```

## License

MIT License - see LICENSE file for details.

## Contributing

Contributions welcome! Please open an issue or pull request.

## Related Packages

- [SpectralModels.jl](../SpectralModels): Underlying spectral fitting framework
- [DmspElfinConjunction.jl](../../): Full analysis pipeline
