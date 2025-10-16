"""
Model evaluation functions for computing fluxes.
"""

"""
    evaluate_flux(model::EmpiricalFluxModel, energy::Real; mlat, mlt, ae, stat=:median)

Evaluate differential flux at a single energy.

# Arguments
- `model`: The empirical flux model
- `energy`: Energy in keV
- `mlat`: Magnetic latitude in degrees
- `mlt`: Magnetic local time in hours (0-24)
- `ae`: AE index in nT
- `stat`: Statistic to use (`:median`, `:q25`, `:q75`)

# Returns
Differential flux in cm⁻² s⁻¹ sr⁻¹ keV⁻¹

# Example
```julia
model = load_flux_model()
flux = evaluate_flux(model, 10.0; mlat=65.0, mlt=6.0, ae=150.0)
```
"""
function evaluate_flux(model::EmpiricalFluxModel, energy::Real;
                       mlat::Real, mlt::Real, ae::Real, stat::Symbol=:median)
    params = interpolate_parameters(model, mlat, mlt, ae; stat=stat)
    spectral_model = to_spectral_model(params)
    return spectral_model(energy)
end

"""
    flux_spectrum(model::EmpiricalFluxModel, energies::AbstractVector;
                  mlat, mlt, ae, stat=:median)

Compute flux spectrum over an energy range.

# Arguments
- `model`: The empirical flux model
- `energies`: Vector of energies in keV
- `mlat`: Magnetic latitude in degrees
- `mlt`: Magnetic local time in hours (0-24)
- `ae`: AE index in nT
- `stat`: Statistic to use (`:median`, `:q25`, `:q75`)

# Returns
Vector of differential flux values in cm⁻² s⁻¹ sr⁻¹ keV⁻¹

# Example
```julia
model = load_flux_model()
E = 10 .^ range(log10(0.03), log10(1000), length=100)
spectrum = flux_spectrum(model, E; mlat=65.0, mlt=6.0, ae=150.0)
```
"""
function flux_spectrum(model::EmpiricalFluxModel, energies::AbstractVector;
                       mlat::Real, mlt::Real, ae::Real, stat::Symbol=:median)
    params = interpolate_parameters(model, mlat, mlt, ae; stat=stat)
    spectral_model = to_spectral_model(params)
    return spectral_model.(energies)
end

"""
    number_flux(model::EmpiricalFluxModel, E_min::Real, E_max::Real;
                mlat, mlt, ae, stat=:median)

Compute integrated number flux over an energy range.

# Arguments
- `model`: The empirical flux model
- `E_min`: Minimum energy in keV
- `E_max`: Maximum energy in keV
- `mlat`: Magnetic latitude in degrees
- `mlt`: Magnetic local time in hours (0-24)
- `ae`: AE index in nT
- `stat`: Statistic to use (`:median`, `:q25`, `:q75`)

# Returns
Number flux in cm⁻² s⁻¹ sr⁻¹

# Example
```julia
model = load_flux_model()
J = number_flux(model, 0.03, 100.0; mlat=65.0, mlt=6.0, ae=150.0)
```
"""
function number_flux(model::EmpiricalFluxModel, E_min::Real, E_max::Real;
                     mlat::Real, mlt::Real, ae::Real, stat::Symbol=:median)
    params = interpolate_parameters(model, mlat, mlt, ae; stat=stat)
    spectral_model = to_spectral_model(params)
    return n_flux(spectral_model, E_min, E_max)
end

"""
    energy_flux(model::EmpiricalFluxModel, E_min::Real, E_max::Real;
                mlat, mlt, ae, stat=:median)

Compute integrated energy flux over an energy range.

# Arguments
- `model`: The empirical flux model
- `E_min`: Minimum energy in keV
- `E_max`: Maximum energy in keV
- `mlat`: Magnetic latitude in degrees
- `mlt`: Magnetic local time in hours (0-24)
- `ae`: AE index in nT
- `stat`: Statistic to use (`:median`, `:q25`, `:q75`)

# Returns
Energy flux in keV cm⁻² s⁻¹ sr⁻¹

# Example
```julia
model = load_flux_model()
JE = energy_flux(model, 0.03, 100.0; mlat=65.0, mlt=6.0, ae=150.0)
```
"""
function energy_flux(model::EmpiricalFluxModel, E_min::Real, E_max::Real;
                     mlat::Real, mlt::Real, ae::Real, stat::Symbol=:median)
    params = interpolate_parameters(model, mlat, mlt, ae; stat=stat)
    spectral_model = to_spectral_model(params)
    return e_flux(spectral_model, E_min, E_max)
end

"""
    flux_components(model::EmpiricalFluxModel, energy::Real;
                    mlat, mlt, ae, stat=:median)

Get individual spectral components at a single energy.

# Returns
Named tuple `(low_energy, high_energy, total)` with flux contributions from:
- `low_energy`: ExpPow model (thermal/core population)
- `high_energy`: Kappa model (suprathermal tail)
- `total`: Combined flux

# Example
```julia
model = load_flux_model()
components = flux_components(model, 50.0; mlat=65.0, mlt=6.0, ae=150.0)
println("Low E: ", components.low_energy)
println("High E: ", components.high_energy)
println("Total: ", components.total)
```
"""
function flux_components(model::EmpiricalFluxModel, energy::Real;
                         mlat::Real, mlt::Real, ae::Real, stat::Symbol=:median)
    params = interpolate_parameters(model, mlat, mlt, ae; stat=stat)
    spectral_model = to_spectral_model(params)

    low_energy = spectral_model.model1(energy)
    high_energy = spectral_model.model2(energy)
    total = low_energy + high_energy

    return (low_energy=low_energy, high_energy=high_energy, total=total)
end
