using JLD2, HDF5

function prepare_time_selection(timestamps, timerange)
    if isnothing(timerange) || !issorted(timestamps)
        time = unix2datetime.(timestamps)
        idxs = nothing
    else
        idx1 = searchsortedfirst(timestamps, datetime2unix(DateTime(timerange[1])))
        idx2 = searchsortedlast(timestamps, datetime2unix(DateTime(timerange[2])))
        idxs = idx1:idx2
        time = @views unix2datetime.(timestamps[idxs])
    end
    return time, idxs
end


function read_and_select(dset::HDF5.Dataset, idxs)
    if HDF5.iscontiguous(dset)
        fdata = HDF5.readmmap(dset)
        return Array(selectdim(fdata, 1, idxs))
    else
        fdata = HDF5.read(dset)
        return selectdim(fdata, 1, idxs)
    end
end

# https://github.com/JuliaIO/JLD2.jl/issues/648
function read_and_select(dset::JLD2.Dataset, idxs)
    fdata = JLD2.readmmap(dset)
    return Array(selectdim(fdata, 1, idxs))
end


function check_hdf_file(path)
    !isfile(path) && error("File not found: $path")
    return jldopen(path) do fid
        # Check if the file has the expected structure
        @info "1D Parameters" fid["Data/Array Layout/1D Parameters"]
        @info "2D Parameters" fid["Data/Array Layout/2D Parameters"]
    end
end

function _read2dimarray(param, params_1d, params_2d, meta_data_params, idxs, all_dims)
    if haskey(params_1d, param)
        dset = params_1d[param]
        dims = all_dims[1]
    elseif haskey(params_2d, param)
        # dset = JLD2.get_dataset(jld_f, "Data/Array Layout/2D Parameters/$param")
        dset = params_2d[param]
        dims = all_dims
    end
    data = isnothing(idxs) ? read(dset) : read_and_select(dset, idxs)
    meta_idx = findfirst(p -> p.mnemonic == uppercase(param), meta_data_params)
    metadata = Dict(pairs(meta_data_params[meta_idx]))
    return DimArray(data, dims; metadata, name = param)
end

# JLD2 is faster and more memory efficient than HDF5 for reading data but not feature complete
function read2dimarray(path, param, timerange = nothing)
    jld_f = jldopen(path)
    fid = h5open(path)

    params_1d = fid["Data/Array Layout/1D Parameters"]
    params_2d = fid["Data/Array Layout/2D Parameters"]
    meta_data_params = read(fid["Metadata/Data Parameters"])

    timestamps = jld_f["Data/Array Layout/timestamps"]
    ch_energy = jld_f["Data/Array Layout/ch_energy"]
    time, idxs = prepare_time_selection(timestamps, timerange)
    dims = (Ti(time), Y(ch_energy))

    da = _read2dimarray(param, params_1d, params_2d, meta_data_params, idxs, dims)
    close(fid)
    close(jld_f)
    return da
end

function read2dimstack(path, params = nothing, timerange = nothing)
    return h5open(path) do fid
        timestamps = read(fid["Data/Array Layout/timestamps"])
        ch_energy = read(fid["Data/Array Layout/ch_energy"])
        params_1d = fid["Data/Array Layout/1D Parameters"]
        params_2d = fid["Data/Array Layout/2D Parameters"]
        meta_data_params = read(fid["Metadata/Data Parameters"])
        params = something(params, setdiff(keys(params_1d), ("Data Parameters",)))

        time, idxs = prepare_time_selection(timestamps, timerange)
        dims = (Ti(time), Y(ch_energy))
        das = map(String.(params)) do param
            _read2dimarray(param, params_1d, params_2d, meta_data_params, idxs, dims)
        end
        metadata = Dict(read(fid["Metadata/Experiment Parameters"]))
        DimStack(das; metadata)
    end
end
