using Tullio, LoopVectorization
using Bumper

include("spectra.jl")


"""
EPD (Energetic Particle Detector) processing functions.

The EPD measures energetic particles with 16 energy channels from ~50 keV to ~5.8 MeV.
Processing includes calibration, flux calculations, and spectral analysis.
"""

# Energy bins for EPD (16 channels, log-spaced from ~50 keV to ~5.8 MeV)
const EPD_ENERGY_BINS = [
    56.28, 75.32, 100.84, 134.99, 180.8, 242.16, 324.42, 434.51,
    582.12, 779.85, 1044.37, 1399.25, 1874.41, 2511.41, 3364.1, 4507.46,
]

unwrap(x) = parent(x)

"""
    epd(trange, probe; type = "nflux", datatype = "pef", fullspin = false)

Load ELFIN EPD CDF data and process ELFIN EPD data to extract flux products (para, anti, precipitating).
"""
function epd(trange, probe; type = "nflux", datatype = "pef", fullspin = false, Espectra = (;), no_update = false, kw...)
    cdf = load(trange, probe; no_update)
    res = fullspin ? :fs : :hs
    base_var = "el$(probe)_$(datatype)_$(res)"
    spec_tvar = "$(base_var)_Epat_$(type)"
    spec_data = cdf[spec_tvar]
    # Get time, pitch angles, and loss cone data
    # time_var = pyconvert(String, cdf.data[spec_tvar].attributes["DEPEND_0"].value)
    # time = cdf[time_var]
    pitch_angles = CDF.dim(spec_data, 2)
    loss_cone = cdf["$(base_var)_LCdeg"]

    S = copy(spec_data)
    PA = copy(pitch_angles)
    sort_flux_by_pitch_angle!(S, PA)
    return merge(
        epd_l2_Espectra(S, PA, loss_cone; fullspin, Espectra...),
        epd_l2_PAspectra(S; kw...)
    )
end

"""
    epd_spectra(flux, pitch_angles, loss_cone; flux_type="nflux", res="hs", para_tol=22.25, anti_tol=22.25)

Process EPD L2 CDF data to create energy spectra for different flux directions from 3D pitch angle data.

# Returns
- DimStack containing omni, para, anti, and prec flux spectra

# Example
```julia
cdf_data = pycdfpp.load(file_path)
spectra = epd_spectra(cdf_data, "a"; flux_type="nflux")
```
"""
function epd_l2_Espectra(flux, pitch_angles, loss_cone; fullspin = false, kw...)
    # Determine nspinsectors based on resolution
    res = fullspin ? :fs : :hs
    n_pa = size(pitch_angles, 2)
    nspinsectors = res == :hs ? (n_pa - 2) * 2 : n_pa - 2

    # Calculate tolerances
    FOVo2 = 11.0  # Field of View divided by 2 (deg)
    half_sector_width = 180 / nspinsectors
    para_tol = FOVo2 + half_sector_width
    perp_tol = -FOVo2
    return Espectra(flux, pitch_angles, loss_cone; para_tol, perp_tol, half_sector_width, kw...)
end

function _sort_by_perm!(x, cache, perm)
    for k in eachindex(x, cache)
        cache[k] = x[perm[k]]
    end
    return x .= cache
end

"""
    sort_flux_by_pitch_angle(flux, pitch_angle)

Sort flux by pitch angle into ascending order. Set `sort_on_reverse` to true to sort only if the pitch angle is in reverse order.
"""
function sort_flux_by_pitch_angle!(flux, pitch_angle; sort_on_reverse = true)
    pa_dim = 2
    n_pa = size(pitch_angle, pa_dim)
    N_flux = ndims(flux)
    N_pa = ndims(pitch_angle)
    @assert size(flux, pa_dim) == n_pa
    @assert N_flux ∈ (2, 3)
    @assert N_pa == 2

    # Sort each time step
    perm = zeros(Int, n_pa)
    T = promote_type(eltype(flux), eltype(pitch_angle))
    cache = zeros(T, n_pa)
    @inbounds for (flux_slice, pa_slice) in zip(eachslice(flux, dims = 1), eachslice(pitch_angle, dims = 1))
        sortperm!(perm, pa_slice)
        sort_on_reverse && perm != n_pa:-1:1 && continue
        _sort_by_perm!(pa_slice, cache, perm)
        for f in eachcol(flux_slice)
            _sort_by_perm!(f, cache, perm)
        end
    end
    return flux, pitch_angle
end

"""
    epd_pitch_angle_spectra(S; energybins=nothing, energies=nothing)

Process EPD L2 CDF data to create pitch angle spectra following Python epd_l2_PAspectra logic.

# Arguments
- `S`: 3D spectral data (time × pitch_angle × energy)
- `energybins`: List of tuples specifying energy channel ranges, e.g., [(0,2), (3,5), (6,8), (9,15)]
- `energies`: List of tuples specifying energy ranges in keV, e.g., [(50,160), (160,345), (345,900), (900,7000)]

# Returns
- Named tuple containing pitch angle spectra for each energy channel
"""
function epd_l2_PAspectra(S; energybins = nothing, energies = nothing)
    # Energy bin boundaries (constant from Python code)
    EMINS = [
        50.000008, 79.999962, 120.00005, 159.99998, 210.00015, 269.99973,
        345.00043, 429.99945, 630.00061, 899.9989, 1300.0013, 1799.9985,
        2500.0022, 3349.999, 4150.0034, 5800.0,
    ]

    EMAXS = [
        79.999962, 120.00005, 159.99998, 210.00015, 269.99973, 345.00043,
        429.99945, 630.00061, 899.9989, 1300.0013, 1799.9985, 2500.0022,
        3349.999, 4150.0034, 5800.0, 7200.0,
    ]
    dE = EMAXS .- EMINS

    # Determine energy channel ranges
    if !isnothing(energybins)
        min_channels = [eb[1] for eb in energybins]
        max_channels = [eb[2] for eb in energybins]
        !isnothing(energies) && @warn "Both energies and energybins are set, 'energybins' takes precedence!"
    elseif !isnothing(energies)
        min_channels = [findfirst(e -> e >= en[1], EPD_ENERGY_BINS) for en in energies]
        max_channels = [findlast(e -> e <= en[2], EPD_ENERGY_BINS) for en in energies]
    else
        # Default energy channels: [(1,3), (4,6), (7,9), (10,16)]
        min_channels = (1, 4, 7, 10)
        max_channels = (3, 6, 9, 16)
    end
    keys = Tuple(Symbol("ch$(i - 1)") for i in 1:length(min_channels))
    values = PAspectra(S, dE, min_channels, max_channels)
    return NamedTuple{keys}(values)
end
