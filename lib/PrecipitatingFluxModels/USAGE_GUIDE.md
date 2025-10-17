### PrecipitatingFluxModels.jl Usage Guide

Complete guide for creating and using empirical flux models.

## Table of Contents

1. [For Model Creators](#for-model-creators)
2. [For End Users](#for-end-users)
3. [Advanced Usage](#advanced-usage)
4. [API Reference](#api-reference)

---

## For Model Creators

### Step 1: Prepare Your Statistical Analysis

Your statistical analysis must produce a DataFrame with the following structure:

```julia
using DataFrames

# Example structure
df = DataFrame(
    # Spatial bins
    mlat_bin = Float64[],      # Magnetic latitude (degrees)
    mlt_bin = Int[],            # Magnetic local time (hours)
    ae_bin = String[],          # AE index bin (e.g., "[0, 100)")
    n_samples = Int[],          # Number of observations in bin

    # For each parameter X, include _median, _q25, _q75:
    # κ, log_A1, E_c1, γ, log_A2, E_c2, Emin, J1, J2, JE1, JE2
    κ_median = Float64[], κ_q25 = Float64[], κ_q75 = Float64[],
    # ... repeat for all parameters
)
```

### Step 2: Create the Model

Use the provided script template:

```julia
using PrecipitatingFluxModels
using DataFrames, DataFramesMeta

# Load your analysis results
include("path/to/your/analysis.jl")  # Should create `sdf` DataFrame

# Prepare parameters (group by spatial bins and compute statistics)
params_df = compute_spatial_statistics(sdf)

# Create model
model = create_model_from_stats(params_df;
    mlat_bins=50:1:80,                              # Your MLat bins
    mlt_bins=0:2:22,                                # Your MLT bins
    ae_bins=["[0, 100)", "[100, 300)", "[300, Inf)"],  # Your AE bins
    energy_range=(0.03, 1000.0),                    # Valid energy range (keV)
    description="DMSP-ELFIN 2020-2022 model",      # Model description
    version="1.0.0"                                  # Version string
)

# Save for distribution
save_flux_model("dmsp_elfin_model_v1.jld2", model)
```

### Step 3: Run the Creation Script

Use the provided script:

```bash
julia scripts/create_flux_model.jl
```

This will:
1. Load your statistical analysis
2. Group observations by spatial bins
3. Compute median and percentiles for all parameters
4. Create the model
5. Save to `lib/PrecipitatingFluxModels/data/default_model.jld2`

---

## For End Users

### Installation

```julia
using Pkg
Pkg.add(url="https://github.com/yourusername/dmsp_elfin")
```

Or for local development:

```julia
Pkg.develop(path="path/to/dmsp_elfin/lib/PrecipitatingFluxModels")
```

### Basic Usage

```julia
using PrecipitatingFluxModels

# Load the model
model = load_flux_model()

# Evaluate flux at a single energy
flux = evaluate_flux(model, 10.0;  # 10 keV
    mlat=65.0,    # Magnetic latitude (degrees)
    mlt=6.0,      # Magnetic local time (hours)
    ae=150.0      # AE index (nT)
)

println("Flux: $flux cm⁻² s⁻¹ sr⁻¹ keV⁻¹")
```

### Computing Energy Spectra

```julia
# Create energy grid (log-spaced)
energies = 10 .^ range(log10(0.03), log10(1000), length=100)

# Compute spectrum
spectrum = flux_spectrum(model, energies;
    mlat=65.0, mlt=6.0, ae=150.0
)

# Plot
using CairoMakie
fig, ax, ln = lines(energies, spectrum;
    axis=(xscale=log10, yscale=log10,
          xlabel="Energy (keV)",
          ylabel="Flux (cm⁻² s⁻¹ sr⁻¹ keV⁻¹)")
)
save("spectrum.png", fig)
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

---

## API Reference

### Model Loading/Saving

```julia
load_flux_model([filename])          # Load model from file or use default
save_flux_model(filename, model)     # Save model to file
create_model_from_stats(df; kwargs)  # Create from statistics DataFrame
```

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

to_spectral_model(params)
    # Convert to TwoStepModel
    # Returns: TwoStepModel

get_ae_bin(ae_value, ae_bins)
    # Find AE bin for value
    # Returns: String
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
3. **Documentation**: Include detailed descriptions of data sources and methodology
4. **Versioning**: Use semantic versioning (MAJOR.MINOR.PATCH)

### For End Users

1. **Interpolation**: Model uses nearest-neighbor interpolation; results most reliable near bin centers
2. **Energy Range**: Respect the model's valid energy range (check `model.energy_range`)
3. **Uncertainty**: Always consider using percentiles for uncertainty estimation
4. **Validation**: Compare model predictions with observations when available

---

## Troubleshooting

### Common Issues

**Model file not found:**
```julia
# Specify explicit path
model = load_flux_model("path/to/model.jld2")
```

**Out of bounds values:**
- Model will use nearest neighbor for out-of-range coordinates
- Extrapolation may not be reliable

**Missing columns error:**
- Ensure all required parameters have _median, _q25, _q75 columns
- Check column names match exactly

---

## Contact

For questions or issues:
- GitHub Issues: [repository URL]
- Email: zijin@ucla.edu
