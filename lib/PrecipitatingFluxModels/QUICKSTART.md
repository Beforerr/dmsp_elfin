# PrecipitatingFluxModels.jl - Quick Start

## For You (Model Creator)

### 1. Generate the Model

After completing your statistical analysis in `notebooks/stats.qmd`:

```bash
cd /Users/zijin/projects/dmsp_elfin
julia scripts/create_flux_model.jl
```

This will:
- Load your `sdf` DataFrame with fitted models
- Group by spatial bins (MLat, MLT, AE)
- Compute median and percentiles
- Create and save the model to `lib/PrecipitatingFluxModels/data/default_model.jld2`

### 2. Test the Package

```bash
julia --project=lib/PrecipitatingFluxModels -e 'using Pkg; Pkg.test()'
```

### 3. Try the Examples

```bash
julia --project=lib/PrecipitatingFluxModels lib/PrecipitatingFluxModels/examples/basic_usage.jl
```

This will generate example plots and demonstrate all features.

---

## For End Users

### Installation

```julia
using Pkg

# Option 1: From GitHub (when published)
Pkg.add(url="https://github.com/yourusername/dmsp_elfin",
        subdir="lib/PrecipitatingFluxModels")

# Option 2: Local development
Pkg.develop(path="/path/to/dmsp_elfin/lib/PrecipitatingFluxModels")
```

### Minimal Example

```julia
using PrecipitatingFluxModels

# Load model
model = load_flux_model()

# Evaluate flux
flux = evaluate_flux(model, 10.0;  # Energy in keV
    mlat=65.0,   # Magnetic latitude (degrees)
    mlt=6.0,     # Magnetic local time (hours, dawn)
    ae=150.0     # AE index (nT, moderate activity)
)

println("Flux at 10 keV: $flux cm⁻² s⁻¹ sr⁻¹ keV⁻¹")
```

### Complete Example with Plot

```julia
using PrecipitatingFluxModels
using CairoMakie

# Load model
model = load_flux_model()

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

## Common Use Cases

### 1. Single Point Evaluation

```julia
flux = evaluate_flux(model, 50.0; mlat=67.0, mlt=8.0, ae=200.0)
```

### 2. Spatial Map

```julia
mlats = 60:1:75
mlts = 0:1:23
flux_map = [evaluate_flux(model, 10.0; mlat=lat, mlt=mlt, ae=150.0)
            for lat in mlats, mlt in mlts]
```

### 3. Activity Scan

```julia
ae_values = 0:50:500
fluxes = [evaluate_flux(model, 10.0; mlat=65.0, mlt=6.0, ae=ae)
          for ae in ae_values]
```

### 4. Integrated Fluxes

```julia
J = number_flux(model, 0.03, 100.0; mlat=65.0, mlt=6.0, ae=150.0)
JE = energy_flux(model, 0.03, 100.0; mlat=65.0, mlt=6.0, ae=150.0)
```

---

## Next Steps

- Read the full [README.md](README.md) for detailed documentation
- Check [USAGE_GUIDE.md](USAGE_GUIDE.md) for advanced patterns
- Run [examples/basic_usage.jl](examples/basic_usage.jl) for comprehensive examples
- See [PACKAGE_SUMMARY.md](PACKAGE_SUMMARY.md) for implementation details

---

## Troubleshooting

### "Model file not found"

Make sure you've generated the model first:
```bash
julia scripts/create_flux_model.jl
```

### "Package not found"

Activate the correct project:
```julia
using Pkg
Pkg.activate("lib/PrecipitatingFluxModels")
```

### Import issues

Make sure SpectralModels is available:
```julia
Pkg.develop(path="/path/to/dmsp_elfin/lib/SpectralModels")
```

---

## Support

For questions or issues:
- GitHub Issues: [repository]
- Email: zijin@ucla.edu
