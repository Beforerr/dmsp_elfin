import Base: read
using URIs

function standardize(x)
    return set(x, Dim{:time} => Ti, Dim{:v_dim} => Y)
end

function apply_date_format(pattern, date)
    return replace(
        pattern,
        "{Y}" => year(date),
        "{M:02d}" => string(month(date); pad = 2),
        "{D:02d}" => string(day(date); pad = 2)
    )
end


# https://github.com/helgee/RemoteFiles.jl/blob/master/src/RemoteFiles.jl
struct RemoteFile
    uri::URI
    path::String
    updates::Symbol
end

function RemoteFile(uri::URI; dir = ".", updates = :never, flatten = false)
    path = local_path(uri, dir; flatten)
    return RemoteFile(uri, path, updates)
end

RemoteFile(uri::String; kwargs...) = RemoteFile(URI(uri); kwargs...)

local_path(f) = f.path

function local_path(uri::URI, dir = "."; flatten::Bool = false)
    fullpath = lstrip(uri.path, '/')
    return if !flatten
        joinpath(dir, fullpath)
    else
        joinpath(dir, basename(fullpath))
    end
end

"""
Set `no_update` to true to only load data from local cache.
"""
function Base.read(f::RemoteFile; no_update::Bool = false)
    path = abspath(local_path(f))
    if no_update
        if isfile(path)
            return read(path)
        end
    else
        any_loc_open = @pyconst pyimport("speasy.core.any_files").any_loc_open
        file = any_loc_open(string(f.uri); cache_remote_files = true)
        return file.read()
    end
end
