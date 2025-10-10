using DataFrames

const FLCS_THRESHOLD = Ref{Float64}(0.5)

@enum Mechanism begin
    FLCS
    WHISTLER
    MIXED
    OTHER
end

"""
    is_nightside(mlt; dusk=20, dawn=4)

Check if MLT is in nightside sector [dusk, 24) ∪ [0, dawn).
"""
is_nightside(mlt; dusk = 20, dawn = 4) = !ismissing(mlt) && (mlt >= dusk || mlt < dawn)

"""
    is_monotone_nonincreasing(vals; tol=1e-6)

Check if values are monotonically non-increasing within tolerance.
"""
is_monotone_nonincreasing(vals; tol = 1.0e-6) =
    length(vals) ≤ 1 || all(≤(tol), diff(vals))

# Replace 0 and Inf with NaN due to noise level
function sanitize(ratio)
    return replace(ratio, 0 => NaN, Inf => NaN)
end

"""
    indicates_flcs(ratio, mlt, Lshell; threshold = FLCS_THRESHOLD)

Check if observation `ratio`: R = |j_para - j_anti| / j_perp` indicates fieldline curvature scattering (FLCS) precipitation.

# Criteria
1. Nightside MLT sector (18-06 MLT)
2. `Lshell` > 5
3. R(maxE) > `threshold`
4. R monotonically decreasing OR staying ≥`threshold` as energy decreases from maxE to energy_floor
"""
function indicates_flcs(ratio, mlt, Lshell = Inf; threshold = 0.5)
    # Check spatial criteria first (fast checks)
    is_nightside(mlt) || return false
    Lshell > 5 || return false
    length(ratio) < 3 && return false
    r_maxE = ratio[end]
    r_maxE < threshold && return false

    # R(100) & R(maxE)>0.75
    return ratio[Energy = Near(100)] > threshold
    # Check monotonic decrease (high to low energy) OR staying above threshold
    # r_descending = reverse(ratio)
    # return is_monotone_nonincreasing(r_descending) || all(≥(threshold), r_descending)
end

# #    63.24554
#    97.97958
#   138.56409
#   183.30309
#   238.11758
#   305.2049
#   385.1623
#   520.48047
#   752.99396
#  1081.6653
#  1529.706
#  2121.3203
#  2893.9602
#  3728.6064
#  4906.1206
#  6500.0
function select_energy(ratio, energy, default = 0)
    range = (energy - 5) .. (energy + 5)
    r = ratio[Energy = range]
    return length(r) > 0 ? r[1] : default
end

"""
    indicates_whistler(ratio; threshold = 0.25)

Check if observation indicates whistler-mode wave precipitation.

# Criteria
1. R(energy_ref) > ratio_threshold
2. R monotonically decreasing OR staying ≥ratio_threshold as energy increases from energy_ref
"""
function indicates_whistler(ratio; threshold = 0.25)
    length(ratio) < 3  && return false
    ratio[1] > threshold || return false
    # Let us try a simple criteria hoping increasing the number
    #  R(63)>1.25*R(200)>1.25*R(300)>1.5*R(500)
    # return select_energy(ratio, 63) > 1.25 * select_energy(ratio, 183) > 1.25 * select_energy(ratio, 305) > 1.5 * select_energy(ratio, 520)
    return is_monotone_nonincreasing(ratio) || ratio[1] > 2 * ratio[Energy = Near(305)]
end

# then 1o MLAT bin can be classified as: FLCS-only (>50% bins satisfy [1] and no bins satisfy [2]), whistler-only (>10% bins satisfy [2] and no bins satisfy [1]), mixed (>25% bins satisfy [1] and >10% bins satisfy [2]), no energetic precipitations (others)

# FLCS -> (more WW) MIXED -> (less FLCS) WHISTLER
function classify(n_total, n_flcs, n_whistler)
    p_flcs = n_flcs / n_total
    p_whistler = n_whistler / n_total
    return if p_flcs >= 0.4 && p_whistler < 0.1
        FLCS
    elseif p_flcs >= 0.25 && p_whistler >= 0.1
        MIXED
    elseif p_flcs < 0.25 && p_whistler >= 0.1
        WHISTLER
    else
        OTHER
    end
end
classify(t) = classify(t...)

function classify_precipitation(ratio, mlt, Lshell = Inf; kw...)
    n_total = size(ratio, 1)
    n_flcs = 0
    n_whistler = 0
    foreach(eachrow(ratio)) do row
        # r = sanitize(row)
        r = row
        r_valid = r[(!isnan).(r)]
        whistler = indicates_whistler(r_valid)
        flcs = !whistler && indicates_flcs(r_valid, mlt)
        n_flcs += flcs
        n_whistler += whistler
    end
    return (n_total, n_flcs, n_whistler)
end

function valid_flux_ratio(flux, trange, threshold)
    prec = tview(flux.prec, trange)
    perp = tview(flux.perp, trange)
    ratio = prec ./ perp
    ratio[prec .< threshold] .= NaN
    ratio[perp .< threshold] .= NaN
    return ratio
end

function classify_precipitation!(df, flux, mlt; trange = :trange, flux_threshold = 100, kw...)
    return @rtransform! df @astable begin
        local_mlt = local_mlt_mean(tview(mlt, $trange))
        ratio = valid_flux_ratio(flux, $trange, flux_threshold)
        :spins_mechanism = classify_precipitation(ratio, local_mlt)
        :mechanism = classify(:spins_mechanism)
    end
end
