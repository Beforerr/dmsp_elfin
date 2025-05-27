module ELFINReader

export read_elfin_data, generate_sample_elfin_data

using CSV
using DataFrames
using Dates
using Downloads
using Random

"""
    read_elfin_data(filename::String)

Read ELFIN data from a CSV file.
Expected columns:
- timestamp: Time in ISO format or unix timestamp
- mlat: Magnetic Latitude in degrees (-90 to 90)
- mlt: Magnetic Local Time in hours (0-24)
"""
function read_elfin_data(filename::String)
    println("Reading ELFIN data from $filename")
    df = CSV.read(filename, DataFrame)
    
    # Ensure timestamp column is datetime
    if "timestamp" in names(df)
        if eltype(df.timestamp) <: Number
            # Convert unix timestamp to datetime
            df.timestamp = unix2datetime.(df.timestamp)
        elseif eltype(df.timestamp) <: AbstractString
            # Convert string timestamp to datetime
            df.timestamp = DateTime.(df.timestamp)
        end
    else
        error("Data must contain a 'timestamp' column")
    end
    
    # Check for required columns
    for col in ["mlt", "mlat"]
        if !(col in names(df))
            error("Data must contain a '$col' column")
        end
    end
    
    return df
end

"""
    generate_sample_elfin_data(dmsp_df::DataFrame, output_file::String)

Generate sample ELFIN data based on DMSP data with some offset to create conjunctions.
"""
function generate_sample_elfin_data(dmsp_df::DataFrame, output_file::String)
    println("Generating sample ELFIN data based on DMSP data")
    
    # Create a copy of DMSP data with some offset to create conjunctions
    elfin_df = copy(dmsp_df)
    
    # Add random offsets to create some conjunctions
    n = nrow(elfin_df)
    
    # Add slight offset to MLT (with some periods of close conjunction)
    Random.seed!(42)  # For reproducibility
    
    # Create MLT offsets that will sometimes be close to DMSP (conjunction)
    # and sometimes be far from DMSP (no conjunction)
    t = range(0, 4π, length=n)
    mlt_offset = 0.5 * sin.(t) .+ 0.2 * sin.(2.5 .* t) .+ 0.1 * randn(n)
    elfin_df.mlt = mod.(elfin_df.mlt .+ mlt_offset, 24)
    
    # Add slight offset to MLAT (with some periods of close conjunction)
    mlat_offset = 2.0 * sin.(range(0, 6π, length=n)) .+ 0.5 * randn(n)
    elfin_df.mlat = elfin_df.mlat .+ mlat_offset
    
    # Save to CSV
    CSV.write(output_file, elfin_df)
    println("Saved sample ELFIN data to $output_file")
    
    return elfin_df
end

end # module
