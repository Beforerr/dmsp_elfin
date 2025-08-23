using Reexport
@reexport using GeoAACGM
import GeoAACGM: geo2aacgm, geod2aacgm
using SPEDAS: times

export gei2aacgm

function GeoAACGM.geo2aacgm(x)
    ts = times(x)
    aacgm = GeoAACGM.geo2aacgm.(eachslice(parent(x), dims = 2)..., ts)
    tdim = Ti(ts)
    mlat = DimArray(getindex.(aacgm, 1), tdim)
    mlon = DimArray(getindex.(aacgm, 2), tdim)
    return DimStack((; mlat, mlon))
end

gei2aacgm(x) = geo2aacgm(gei2geo(x))

function GeoAACGM.geod2aacgm(x::DimStack)
    times = parent(dims(x, Ti))
    mlat, mlon, r = similar.(values(x))
    for i in eachindex(x)
        mlat[i], mlon[i], r[i] = geod2aacgm(x[i]..., times[i])
    end
    return DimStack((; mlat, mlon, r))
end
