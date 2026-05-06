using SpectralModels
using Chairmarks
import SpectralModels.NonlinearSolve as NonlinearSolve

# Data from event id=16, probe=a, 2022-07-02T05:00:00
const E_dmsp = [
    0.03, 0.044, 0.065, 0.095, 0.139, 0.20400000000000001, 0.3, 0.44, 0.646,
    0.9490000000000001, 1.3920000000000001, 2.04, 3.0, 4.4, 6.46,
    9.450000000000001, 13.9, 20.400000000000002, 30.0,
]
const flux_dmsp = [
    3.7792894935752077e9, 8.98876404494382e9, 1.0680447889750216e10,
    1.208533727956591e10, 6.301824212271974e9, 4.728132387706855e9,
    5.39452495974235e9, 4.173716864072194e9, 4.8469911686177845e9,
    3.790310255234601e9, 2.869471413160733e9, 2.70735524256651e9,
    1.9778941244909828e9, 1.0941475826972008e9, 3.1200254695956707e8,
    1.4556040756914118e8, 7.382798080472499e6, 0.0, 4.257130693912302e6,
]

const E_elx = [
    63.24554, 97.97958, 138.56409, 183.30309, 238.11758, 305.2049,
    385.1623, 520.48047, 752.99396, 1081.6653, 1529.706, 2121.3203,
    2893.9602, 3728.6064, 4906.1206, 6500.0,
]
const flux_elx = [
    16003.93, 2658.736, 1364.682, 1111.4905, 93.81907, 75.05525,
    0.0, 0.0, 47.725285, 7.0531034, 12.996939, 0.0, 6.1931925,
    7.255754, 0.8348091, 8.490667,
]

const THRESH = 150.0
const ModelType = TwoStepModel{PowerLawExpCutoff2, TransformKappaDistribution}

E = vcat(E_dmsp, E_elx)
flux = vcat(flux_dmsp, flux_elx)
valid = @. !isnan(flux) & (flux >= THRESH)
E_clean = E[valid]
flux_clean = flux[valid]

result, score = SpectralModels.fit(ModelType, flux_clean, E_clean)
result2, score2 = SpectralModels.fit(ModelType, flux_clean, E_clean; refine = true)
@info result score
@info result2 score2
@info @b SpectralModels.fit(ModelType, flux_clean, E_clean)
@info @b SpectralModels.fit(ModelType, flux_clean, E_clean; refine = true)
