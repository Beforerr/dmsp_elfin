import GeoAACGM
import GeoAACGM: geo2aacgm, geod2aacgm, geod2geoc, geod2geo
using SPEDAS: times

export gei2aacgm, geo2aacgm, geod2aacgm
export geod2geoc, geod2geo

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

function GeoAACGM.geod2geoc(x::DimStack)
    lat, lon, height = similar.(values(x))
    for i in eachindex(x)
        lat[i], lon[i], height[i] = geod2geoc(x[i]...)
    end
    return DimStack((; lat, lon, height))
end


function GeoAACGM.geod2geo(geod::DimStack)
    x, y, z = similar.(values(geod))
    for i in eachindex(geod)
        x[i], y[i], z[i] = geod2geo(geod[i]...)
    end
    return DimStack((; x, y, z))
end
