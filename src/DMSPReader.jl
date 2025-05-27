module DMSPReader

export read_dmsp_ssj_file, download_dmsp_data, read_dmsp_ssj_records

using DataFrames
using Dates
using Downloads
using GZip

"""
    download_dmsp_data(url::String, output_file::String)

Download DMSP data from the specified URL and save it to the output file.
"""
function download_dmsp_data(url::String, output_file::String)
    println("Downloading DMSP data from $url")
    Downloads.download(url, output_file)
    println("Downloaded to $output_file")
    return output_file
end


"""
    DMSPRecord

Structure representing a single DMSP data record with processed values.
"""
struct DMSPRecord
    timestamp::DateTime
    lat::Float64              # Geodetic latitude in degrees
    lon::Float64              # Geographic longitude in degrees
    altitude::Float64         # Altitude in nautical miles
    mlat::Float64             # Corrected geomagnetic latitude in degrees
    mlon::Float64             # Corrected geomagnetic longitude in degrees
    mlt::Float64              # Magnetic local time in hours
    electron_flux::Vector{Float64}  # Electron flux for 20 channels
    ion_flux::Vector{Float64}      # Ion flux for 20 channels
    electron_energy::Vector{Float64}  # Electron energy levels in eV
    ion_energy::Vector{Float64}      # Ion energy levels in eV
end


# Structure to hold data for a single second
struct SecondData
    hour::UInt16               # Hour of day for this second
    minute::UInt16             # Minute of hour for this second
    second::UInt16             # Second of minute for this second
    electron_flux::NTuple{20,UInt16}  # Electron flux data for 20 channels
    ion_flux::NTuple{20,UInt16}   # Ion flux data for 20 channels
end


"""
    SSJRecord

Raw binary structure of a DMSP SSJ record based on the new file format description.

The SSJ data files are binary files with the data stored as a series of 16-bit unsigned integers
written using big endian encoding.

File format:
- Word 1: Day of year, 1 to 366, days (2 bytes)
- Word 2: Hour of day, 0 to 23, hours (2 bytes)
- Word 3: Minute of hour, 0 to 59, minutes (2 bytes)
- Word 4: Second of minute, 0 to 59, seconds (2 bytes)
- Word 5: Integer year, 1987 to 2049, years, conversion: i - 50 (2 bytes)
- Word 6: Geodetic latitude, -90.0 to 90.0, degrees, conversion: float(i-900)/10.0 if i > 1800 then float(i - 4995)/10.0 (2 bytes)
- ... (remaining fields as per the format description)
"""
struct SSJRecord
    # Header information (words 1-15)
    day_of_year::UInt16       # Word 1: Day of year, 1-366
    hour::UInt16              # Word 2: Hour of day, 0-23
    minute::UInt16            # Word 3: Minute of hour, 0-59
    second::UInt16            # Word 4: Second of minute, 0-59
    year_code::UInt16         # Word 5: Integer year (need to subtract 50 to get actual year)
    geodetic_lat_code::UInt16 # Word 6: Geodetic latitude code
    geo_lon_code::UInt16      # Word 7: Geographic longitude code
    altitude::UInt16          # Word 8: Altitude in nautical miles
    geo_lat_110km_code::UInt16 # Word 9: Geographic latitude at 110 km altitude code
    geo_lon_110km_code::UInt16 # Word 10: Geographic longitude at 110 km altitude code
    cgm_lat_code::UInt16      # Word 11: Corrected geomagnetic latitude code
    cgm_lon_code::UInt16      # Word 12: Corrected geomagnetic longitude code
    mlt_hour::UInt16          # Word 13: Hour of magnetic local time
    mlt_minute::UInt16        # Word 14: Minute of magnetic local time
    mlt_second::UInt16        # Word 15: Second of magnetic local time

    # All seconds of data (60 seconds per minute)
    seconds_data::Vector{SecondData}  # Data for all seconds in this minute

    # Flag for millisecond precision
    ms_flag::UInt16           # Word 2596: Flag for millisecond precision
end

nread(io, type) = ntoh(read(io, type))

# Function to read a single second of data
function read_second_data(io::IO)
    # Read time information for this second (words 16-18)
    hour = nread(io, UInt16)         # Word 16: Hour of day
    minute = nread(io, UInt16)       # Word 17: Minute of hour
    second = nread(io, UInt16)       # Word 18: Second of minute
    # Read electron flux data (words 19-38) and ion flux data (words 39-58)
    electron_flux_array = ntuple(i -> nread(io, UInt16), 20)
    ion_flux_array = ntuple(i -> nread(io, UInt16), 20)
    return SecondData(hour, minute, second, electron_flux_array, ion_flux_array)
end

function read_ssj_record(io::IO)
    # Read header information (words 1-15)
    day_of_year = nread(io, UInt16)       # Word 1: Day of year
    hour = nread(io, UInt16)              # Word 2: Hour of day
    minute = nread(io, UInt16)            # Word 3: Minute of hour
    second = nread(io, UInt16)            # Word 4: Second of minute
    year_code = nread(io, UInt16)         # Word 5: Integer year code
    geodetic_lat_code = nread(io, UInt16) # Word 6: Geodetic latitude code
    geo_lon_code = nread(io, UInt16)      # Word 7: Geographic longitude code
    altitude = nread(io, UInt16)          # Word 8: Altitude in nautical miles
    geo_lat_110km_code = nread(io, UInt16) # Word 9: Geographic latitude at 110 km altitude code
    geo_lon_110km_code = nread(io, UInt16) # Word 10: Geographic longitude at 110 km altitude code
    cgm_lat_code = nread(io, UInt16)      # Word 11: Corrected geomagnetic latitude code
    cgm_lon_code = nread(io, UInt16)      # Word 12: Corrected geomagnetic longitude code
    mlt_hour = nread(io, UInt16)          # Word 13: Hour of magnetic local time
    mlt_minute = nread(io, UInt16)        # Word 14: Minute of magnetic local time
    mlt_second = nread(io, UInt16)        # Word 15: Second of magnetic local time
    # Read all 60 seconds of data
    # Each second of data is 43 words (words 16-58)
    seconds_data = map(i -> read_second_data(io), 1:60)

    # Read millisecond flag (word 2596)
    ms_flag = nread(io, UInt16)           # Word 2596: Flag for millisecond precision
    # Skip the rest of the data for this minute (words 2597-2640)
    skip(io, 2 * (2640 - 2596))  # Skip to the end of this minute's data, 2 bytes per word

    return SSJRecord(
        day_of_year, hour, minute, second, year_code, geodetic_lat_code,
        geo_lon_code, altitude, geo_lat_110km_code, geo_lon_110km_code,
        cgm_lat_code, cgm_lon_code, mlt_hour, mlt_minute, mlt_second,
        seconds_data, ms_flag
    )
end

"""
    read_dmsp_ssj_records(filename)

Read all SSJ records from a file and return them as a vector of SSJRecord.
"""
function read_dmsp_ssj_records(filename)
    # Open the file (handle both .gz and uncompressed)
    if endswith(filename, ".gz")
        io = GZip.open(filename)
    else
        io = open(filename, "r")
    end
    # Calculate number of records based on the format
    # Each minute of data is 2640 words (5280 bytes)
    # The file size should be a multiple of 5280 bytes
    file_size = filesize(filename)
    bytes_per_minute = 2640 * 2  # 2640 words * 2 bytes per word
    num_minutes = div(file_size, bytes_per_minute)
    @info "File contains data for approximately $num_minutes minutes"
    # Read all minutes of data
    records = map(i -> read_ssj_record(io), 1:num_minutes)
    close(io)
    return records
end

"""
    convert_to_dmsp_record(record::SSJRecord)

Convert an SSJRecord to a DMSPRecord with physical units using the new format conversion rules.
"""
function convert_to_dmsp_record(record::SSJRecord, second_index::Int=1)
    # Get the data for the specified second
    second_data = record.seconds_data[second_index]

    # Convert year code to actual year (i - 50)
    year = record.year_code - 50

    # Create timestamp using the header date and the time from the specific second
    date = Dates.DateTime(year, 1, 1) + Dates.Day(record.day_of_year - 1) +
           Dates.Hour(second_data.hour) + Dates.Minute(second_data.minute) + Dates.Second(second_data.second)

    # Convert geodetic latitude code to actual latitude
    lat = if record.geodetic_lat_code > 1800
        (record.geodetic_lat_code - 4995) / 10.0
    else
        (record.geodetic_lat_code - 900) / 10.0
    end

    # Convert geographic longitude code to actual longitude
    lon = record.geo_lon_code / 10.0

    # Convert altitude (already in nautical miles)
    altitude = Float64(record.altitude)

    # Convert magnetic coordinates
    mlat = if record.cgm_lat_code > 1800
        (record.cgm_lat_code - 4995) / 10.0
    else
        (record.cgm_lat_code - 900) / 10.0
    end

    mlon = record.cgm_lon_code / 10.0

    # Calculate MLT (magnetic local time)
    mlt = record.mlt_hour + record.mlt_minute / 60.0 + record.mlt_second / 3600.0

    # Define energy levels for each channel (in eV)
    electron_energy = [
        9450.0, 13900.0, 20400.0, 30000.0,  # Channels 4-1
        2040.0, 3000.0, 4400.0, 6460.0,      # Channels 8-5
        646.0, 949.0, 949.0, 1392.0,         # Channels 12-9
        139.0, 204.0, 300.0, 440.0,          # Channels 16-13
        30.0, 44.0, 65.0, 95.0               # Channels 20-17
    ]

    ion_energy = [
        9450.0, 13900.0, 20400.0, 30000.0,  # Channels 4-1
        2040.0, 3000.0, 4400.0, 6460.0,      # Channels 8-5
        646.0, 949.0, 949.0, 1392.0,         # Channels 12-9
        139.0, 204.0, 300.0, 440.0,          # Channels 16-13
        30.0, 44.0, 65.0, 95.0               # Channels 20-17
    ]

    # Convert electron and ion flux to physical values from the specific second
    # Note: The actual conversion might depend on calibration factors
    # For now, we'll just copy the raw values
    e_flux = Float64.([second_data.electron_flux...])  # Convert tuple to vector
    i_flux = Float64.([second_data.ion_flux...])      # Convert tuple to vector

    # Create a DMSPRecord with all the processed data
    return DMSPRecord(
        date, lat, lon, altitude, mlat, mlon, mlt,
        e_flux, i_flux, electron_energy, ion_energy
    )
end

"""
    read_dmsp_ssj_file(filename::String)

Read a DMSP SSJ binary file and return a DataFrame with the data.
"""
function read_dmsp_ssj_file(filename::String; include_all_seconds::Bool=false)
    # Read all records
    records = read_dmsp_ssj_records(filename)

    # Convert to DMSPRecords - either just the first second of each minute or all seconds
    dmsp_records = Vector{DMSPRecord}()

    # Include all 60 seconds from each minute
    for record in records
        for second_index in 1:60
            push!(dmsp_records, convert_to_dmsp_record(record, second_index))
        end
    end

    # Extract data for DataFrame
    timestamps = [record.timestamp for record in dmsp_records]
    lats = [record.lat for record in dmsp_records]
    lons = [record.lon for record in dmsp_records]
    altitudes = [record.altitude for record in dmsp_records]
    mlats = [record.mlat for record in dmsp_records]
    mlons = [record.mlon for record in dmsp_records]
    mlts = [record.mlt for record in dmsp_records]

    # Create DataFrame with basic fields
    df = DataFrame(
        timestamp=timestamps,
        lat=lats,
        lon=lons,
        altitude=altitudes,
        mlat=mlats,
        mlon=mlons,
        mlt=mlts
    )

    # Add electron and ion flux columns
    if !isempty(dmsp_records) && !isempty(dmsp_records[1].electron_flux)
        # Get the energy levels for column names
        electron_energies = dmsp_records[1].electron_energy
        ion_energies = dmsp_records[1].ion_energy

        # Add electron flux columns
        for i in 1:length(dmsp_records[1].electron_flux)
            energy = electron_energies[i]
            df[!, Symbol("electron_flux_$(Int(energy))eV")] = [record.electron_flux[i] for record in dmsp_records]
        end

        # Add ion flux columns
        for i in 1:length(dmsp_records[1].ion_flux)
            energy = ion_energies[i]
            df[!, Symbol("ion_flux_$(Int(energy))eV")] = [record.ion_flux[i] for record in dmsp_records]
        end
    end

    return df
end

end