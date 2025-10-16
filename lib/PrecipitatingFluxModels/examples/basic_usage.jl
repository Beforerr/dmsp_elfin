using PrecipitatingFluxModels
using CairoMakie

# ============================================================================
# Example 3: Computing an Energy Spectrum
# ============================================================================
println("\n📈 Example 3: Energy Spectrum")
println("-"^60)

# Create log-spaced energy grid
energies = 10 .^ range(log10(0.03), log10(1000), length = 100)
spectrum = flux_spectrum(model, energies; mlat = mlat, mlt = mlt, ae = ae)

println("Computed spectrum over $(length(energies)) energy points")
println("Energy range: $(energies[1]) - $(energies[end]) keV")

# Plot the spectrum
fig = Figure(size = (800, 600))
ax = Axis(
    fig[1, 1];
    xlabel = "Energy (keV)",
    ylabel = "Differential Flux (cm⁻² s⁻¹ sr⁻¹ keV⁻¹)",
    xscale = log10,
    yscale = log10,
    title = "Precipitating Flux Spectrum\nMLat=$(mlat)°, MLT=$(mlt)h, AE=$(ae) nT"
)
lines!(ax, energies, spectrum; linewidth = 2, color = :steelblue)
save("example_spectrum.png", fig)
println("✓ Saved spectrum plot to example_spectrum.png")

# ============================================================================
# Example 4: Uncertainty Quantification
# ============================================================================
println("\n📊 Example 4: Uncertainty Quantification")
println("-"^60)

flux_median = evaluate_flux(model, energy; mlat = mlat, mlt = mlt, ae = ae, stat = :median)
flux_q25 = evaluate_flux(model, energy; mlat = mlat, mlt = mlt, ae = ae, stat = :q25)
flux_q75 = evaluate_flux(model, energy; mlat = mlat, mlt = mlt, ae = ae, stat = :q75)

println("Flux at $(energy) keV with uncertainty:")
println("  Median: $(flux_median)")
println("  25th percentile: $(flux_q25)")
println("  75th percentile: $(flux_q75)")
println("  Range: $(flux_q75 - flux_q25) (IQR)")

# ============================================================================
# Example 5: Integrated Fluxes
# ============================================================================
println("\n∫ Example 5: Integrated Fluxes")
println("-"^60)

E_min, E_max = 0.03, 100.0

J = number_flux(model, E_min, E_max; mlat = mlat, mlt = mlt, ae = ae)
JE = energy_flux(model, E_min, E_max; mlat = mlat, mlt = mlt, ae = ae)

println("Integrated over $(E_min)-$(E_max) keV:")
println("  Number flux (J): $(J) cm⁻² s⁻¹ sr⁻¹")
println("  Energy flux (JE): $(JE) keV cm⁻² s⁻¹ sr⁻¹")

# ============================================================================
# Example 6: Spectral Components
# ============================================================================
println("\n🔬 Example 6: Spectral Components")
println("-"^60)

components = flux_components(model, 50.0; mlat = mlat, mlt = mlt, ae = ae)

println("Flux contributions at 50 keV:")
println("  Low-energy (ExpPow): $(components.low_energy)")
println("  High-energy (Kappa): $(components.high_energy)")
println("  Total: $(components.total)")
println("  Kappa fraction: $(round(components.high_energy / components.total * 100, digits = 1))%")

# ============================================================================
# Example 7: Spatial Mapping
# ============================================================================
println("\n🗺️  Example 7: Spatial Flux Map")
println("-"^60)

mlats = 60:1:75
mlts = 0:1:23
energy_map = 10.0  # keV
ae_map = 150.0  # nT

flux_map = [
    evaluate_flux(model, energy_map; mlat = mlat, mlt = mlt, ae = ae_map)
        for mlat in mlats, mlt in mlts
]

fig2 = Figure(size = (1000, 600))
ax2 = Axis(
    fig2[1, 1];
    xlabel = "MLT (hours)",
    ylabel = "MLat (degrees)",
    title = "$(energy_map) keV Flux Spatial Distribution (AE=$(ae_map) nT)"
)
hm = heatmap!(
    ax2, collect(mlts), collect(mlats), flux_map';
    colorscale = log10,
    colormap = :viridis
)
Colorbar(fig2[1, 2], hm; label = "Flux (cm⁻² s⁻¹ sr⁻¹ keV⁻¹)")
save("example_spatial_map.png", fig2)
println("✓ Saved spatial map to example_spatial_map.png")
println("Flux range: $(minimum(flux_map)) - $(maximum(flux_map))")

# ============================================================================
# Example 8: AE Dependence
# ============================================================================
println("\n⚡ Example 8: Geomagnetic Activity Dependence")
println("-"^60)

ae_values = [50.0, 150.0, 400.0]  # Quiet, moderate, active
energy_ae = 10.0  # keV

println("Flux at $(energy_ae) keV for different activity levels:")
for ae_val in ae_values
    flux_ae = evaluate_flux(model, energy_ae; mlat = mlat, mlt = mlt, ae = ae_val)
    println("  AE=$(ae_val) nT: $(flux_ae)")
end

# ============================================================================
# Example 9: Model Parameters
# ============================================================================
println("\n⚙️  Example 9: Accessing Model Parameters")
println("-"^60)

params = flux_parameters(model; mlat = mlat, mlt = mlt, ae = ae)

println("Model parameters at MLat=$(mlat)°, MLT=$(mlt)h, AE=$(ae) nT:")
println("  κ: $(params.κ)")
println("  E_c1 (ExpPow): $(params.E_c1) keV")
println("  E_c2 (Kappa): $(params.E_c2) keV")
println("  γ: $(params.γ)")
println("  Emin (transition): $(params.Emin) keV")
println("  Number of samples: $(params.n_samples)")

println("\n" * "="^60)
println("✓ All examples completed successfully!")
println("="^60)
