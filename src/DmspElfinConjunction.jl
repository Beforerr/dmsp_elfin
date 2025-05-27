module DmspElfinConjunction

using Dates
using HDF5
using DimensionalData

export read2dimstack

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
