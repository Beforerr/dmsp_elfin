import Base: read
using URIs

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
