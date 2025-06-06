module DmspElfinConjunction

using Dates
using HDF5
using DimensionalData

export read2dimarray, read2dimstack

function read2dimarray(path, param, timerange=nothing)
    h5open(path) do fid
        params_1d = fid["Data/Array Layout/1D Parameters"]
        params_2d = fid["Data/Array Layout/2D Parameters"]
        meta_data_params = read(fid["Metadata/Data Parameters"])
        meta_idx = findfirst(p -> p.mnemonic == uppercase(param), meta_data_params)
        metadata = Dict(pairs(meta_data_params[meta_idx]))

        if haskey(params_1d, param)
            dset = params_1d[param]
        elseif haskey(params_2d, param)
            dset = params_2d[param]
        end
        ndim = HDF5.ndims(dset)

        if ndim == 2
            ch_energy = read(fid["Data/Array Layout/ch_energy"])
            y = Y(ch_energy)
        end

        timestamps = read(fid["Data/Array Layout/timestamps"])
        if isnothing(timerange) || !issorted(timestamps)
            time = Ti(unix2datetime.(timestamps))
            data = read(dset)
            dims = ndim == 1 ? (time,) : (time, y)
        else
            idx1 = findfirst(>=(datetime2unix(timerange[1])), timestamps)
            idx2 = findfirst(>=(datetime2unix(timerange[2])), timestamps)
            time = @views Ti(unix2datetime.(timestamps[idx1:idx2]))
            if HDF5.ismmappable(dset)
                fdata = HDF5.readmmap(dset)
                data = Array(selectdim(fdata, 1, idx1:idx2))
            else
                fdata = read(dset)
                data = selectdim(fdata, 1, idx1:idx2)
            end
            dims = ndim == 1 ? (time,) : (time, y)
        end
        DimArray(data, dims; metadata, name=param)
    end
end

function read2dimstack(path, params=nothing)
    h5open(path) do fid
        timestamps = read(fid["Data/Array Layout/timestamps"])
        time = Ti(unix2datetime.(timestamps))
        ch_energy = read(fid["Data/Array Layout/ch_energy"])
        y = Y(ch_energy)
        params_1d = fid["Data/Array Layout/1D Parameters"]
        params_2d = fid["Data/Array Layout/2D Parameters"]
        meta_data_params = read(fid["Metadata/Data Parameters"])
        params = something(params, setdiff(keys(params_1d), ("Data Parameters",)))

        das = map(String.(params)) do param
            meta_idx = findfirst(p -> p.mnemonic == uppercase(param), meta_data_params)
            metadata = Dict(pairs(meta_data_params[meta_idx]))
            if haskey(params_1d, param)
                data = read(params_1d[param])
                dims = (time,)
            elseif haskey(params_2d, param)
                data = read(params_2d[param])
                # Replace zeros with NaN for 2D parameters
                dims = (time, y)
            end
            DimArray(data, dims; metadata, name=param)
        end
        metadata = Dict(read(fid["Metadata/Experiment Parameters"]))
        DimStack(das; metadata)
    end
end

end
