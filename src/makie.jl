import Makie
import Makie: convert_arguments
using Makie: PointBased

"""
    convert_arguments(::PointBased, energies, model::SpectralModel)

Makie conversion that allows plotting spectral models directly.

# Usage
```julia
E = logrange(10, 1000, 100)
model = PowerLaw(1000.0, 2.5)
lines(E, model)  # Automatically converts to lines(E, model.(E))
```
"""
function Makie.convert_arguments(::PointBased, energies, model::SpectralModel)
    return (energies, model.(energies))
end