# Spectral Models

This page demonstrates the different spectral models available for fitting particle flux data. These models are commonly used in space physics to characterize energy distributions of charged particles.

```@example flux
using DmspElfinConjunction

# Use CairoMakie in CI, GLMakie otherwise
if haskey(ENV, "CI") && ENV["CI"] == "true"
    using CairoMakie
    CairoMakie.activate!()
else
    using GLMakie
    GLMakie.activate!()
end
```

```@example flux
E = logrange(1e-1, 1e3, length=100)
axis = (
    xlabel = "Energy (keV)", 
    ylabel = "Flux", 
    xscale = log10, 
    yscale = log10
)
```

## Power Law

Simple power law model: $f(E) = A · E^{-γ}$

```@example flux
model = PowerLaw(1000.0, 2.5)
lines(E, model; axis, label="PowerLaw")
```

## Power Law with Exponential Cutoff

Power law with exponential cutoff: $f(E) = A · E^{-γ} · \exp(-E/E_c)$

```@example flux
model1 = PowerLawExpCutoff(1.0, 2.5, 10.0)
model2 = PowerLawExpCutoff(1.0, 2.5, 100.0)
model3 = PowerLawExpCutoff(1.0, -1.0, 10.0)

f = lines(E, model1; axis, label="Ec=10 keV")
lines!(E, model2; label="Ec=100 keV") 
lines!(E, model3; label="γ=-1.0")
axislegend()
ylims!(1e-7, 1e3)
f
```

The exponential cutoff energy Ec controls where the spectrum rolls over. A negative γ creates an inverted spectrum that peaks before cutting off.

## Kappa Distribution

Kappa distribution model: $f(E) = A · E · (1 + E/(κ · E_c))^{-κ-1}$

```@example flux
model1 = KappaDistribution(1.0, 2.0, 10.0)
model2 = KappaDistribution(1.0, 5.0, 10.0)
model3 = KappaDistribution(1.0, 10.0, 10.0)
model4 = KappaDistribution(1.0, 10.0, 1.0)

f = lines(E, model1; axis, label="κ=2")
lines!(E, model2; label="κ=5")
lines!(E, model3; label="κ=10")
lines!(E, model4; label="κ=40")
axislegend()
f
```

The kappa parameter controls the suprathermal tail: lower κ values create more pronounced high-energy tails, while higher κ values approach Maxwellian distributions.

### Fitting Kappa Distribution

```@example flux
test = vcat(flux1, flux2)[Y=Where(>=(7))]
test = test[test.>200]
out = DmspElfinConjunction.fit(KappaDistribution, test)
# out = KappaDistribution(init_guess(KappaDistribution, test)...)
f = plot(test; axis)
lines!(DmspElfinConjunction.energies(test), out)
f
```

## Smooth Broken Power Law

Smooth broken power law: $f(E) = A · E^{-γ_1} · (1 + (E/E_b)^m)^{(γ_2 - γ_1)/m}$

```@example flux
model1 = SmoothBrokenPowerlaw(1.0, 1.5, 3.5, 100.0, 1.0)
model2 = SmoothBrokenPowerlaw(1.0, 1.5, 3.5, 200.0, 1.0)
model3 = SmoothBrokenPowerlaw(1.0, 0.5, 4.0, 150.0, 1.0)

f = lines(E, model1; axis, label="Eb=100 keV")
lines!(E, model2; label="Eb=200 keV")
lines!(E, model3; label="γ₁=0.5, γ₂=4.0")
axislegend()
f
```

This model smoothly transitions from a power law with index γ₁ at low energies to γ₂ at high energies, with the transition occurring around the break energy Eb.