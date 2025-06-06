module DmspElfinConjunction

using Dates
using HDF5
using JLD2
using DimensionalData

export read2dimarray, read2dimstack

include("DMSP.jl")

function read_and_select(dset, idxs)
    if HDF5.ismmappable(dset)
        # https://github.com/JuliaIO/JLD2.jl/issues/648
        fdata = HDF5.readmmap(dset)
        data = Array(selectdim(fdata, 1, idxs))
    else
        fdata = read(dset)
        data = selectdim(fdata, 1, idxs)
    end
    return data
end

function prepare_time_selection(timestamps, timerange)
    if isnothing(timerange) || !issorted(timestamps)
        time = unix2datetime.(timestamps)
        idxs = nothing
    else
        idx1 = searchsortedfirst(timestamps, datetime2unix(timerange[1]))
        idx2 = searchsortedfirst(timestamps, datetime2unix(timerange[2]))
        idxs = idx1:idx2
        time = @views unix2datetime.(timestamps[idxs])
    end
    return time, idxs
end

# JLD2 is faster and more memory efficient than HDF5 for reading data but not feature complete
function read2dimarray(path, param, timerange = nothing)
    jld_f = jldopen(path)
    fid = h5open(path)

    params_1d = fid["Data/Array Layout/1D Parameters"]
    params_2d = fid["Data/Array Layout/2D Parameters"]
    meta_data_params = read(fid["Metadata/Data Parameters"])
    meta_idx = findfirst(p -> p.mnemonic == uppercase(param), meta_data_params)
    metadata = Dict(pairs(meta_data_params[meta_idx]))

    timestamps = jld_f["Data/Array Layout/timestamps"]
    time, idxs = prepare_time_selection(timestamps, timerange)

    if haskey(params_1d, param)
        dset = params_1d[param]
        dims = (Ti(time),)
    elseif haskey(params_2d, param)
        ch_energy = jld_f["Data/Array Layout/ch_energy"]
        dset = params_2d[param]
        dims = (Ti(time), Y(ch_energy))
    end
    data = isnothing(idxs) ? read(dset) : read_and_select(dset, idxs)
    da = DimArray(data, dims; metadata, name = param)
    close(fid)
    close(jld_f)
    return da
end

function read2dimstack(path, params = nothing, timerange = nothing)
    return h5open(path) do fid
        timestamps = read(fid["Data/Array Layout/timestamps"])
        ch_energy = read(fid["Data/Array Layout/ch_energy"])
        y = Y(ch_energy)
        params_1d = fid["Data/Array Layout/1D Parameters"]
        params_2d = fid["Data/Array Layout/2D Parameters"]
        meta_data_params = read(fid["Metadata/Data Parameters"])
        params = something(params, setdiff(keys(params_1d), ("Data Parameters",)))

        time, idxs = prepare_time_selection(timestamps, timerange)
        time = Ti(time)

        das = map(String.(params)) do param
            meta_idx = findfirst(p -> p.mnemonic == uppercase(param), meta_data_params)
            metadata = Dict(pairs(meta_data_params[meta_idx]))
            if haskey(params_1d, param)
                dset = params_1d[param]
                dims = (time,)
            elseif haskey(params_2d, param)
                dset = params_2d[param]
                dims = (time, y)
            end
            data = isnothing(idxs) ? read(dset) : read_and_select(dset, idxs)
            DimArray(data, dims; metadata, name = param)
        end
        metadata = Dict(read(fid["Metadata/Experiment Parameters"]))
        DimStack(das; metadata)
    end
end

end
