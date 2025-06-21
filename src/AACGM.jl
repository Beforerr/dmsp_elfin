using Reexport
@reexport using GeoAACGM
import GeoAACGM: geo2aacgm, geod2aacgm
using SPEDAS: times

function GeoAACGM.geo2aacgm(x)
    ts = times(x)
    aacgm = GeoAACGM.geo2aacgm.(eachslice(x, dims=2)..., ts)
    tdim = Ti(ts)
    mlat = DimArray(getindex.(aacgm, 1), tdim)
    mlon = DimArray(getindex.(aacgm, 2), tdim)
    return DimStack((; mlat, mlon))
end

function GeoAACGM.geod2aacgm(x)
    ts = times(x)
    aacgm = GeoAACGM.geod2aacgm.(eachslice(x, dims=2)..., ts)
    tdim = Ti(ts)
    mlat = DimArray(getindex.(aacgm, 1), tdim)
    mlon = DimArray(getindex.(aacgm, 2), tdim)
    r = DimArray(getindex.(aacgm, 3), tdim)
    return DimStack((; mlat, mlon, r))
end

gei2aacgm(x) = geo2aacgm(gei2geo(x))