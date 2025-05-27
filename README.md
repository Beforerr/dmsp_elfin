# DMSP-ELFIN Conjunction Finder

This tool identifies conjunction events between DMSP and ELFIN satellites based on their MLT (Magnetic Local Time) and MLAT (Magnetic Latitude) values.

## What is a Conjunction Event?

A conjunction event is defined when:
- The absolute difference in MLT between DMSP and ELFIN is less than a threshold d₁
- The absolute difference in MLAT between DMSP and ELFIN is less than a threshold d₂
- This condition persists for at least Tₘᵢₙ seconds

## Usage

### Basic Usage

```bash
julia conjunction_finder.jl --dmsp path/to/dmsp_data.gz --elfin path/to/elfin_data.csv
```

### Generate Sample Data

If you don't have real ELFIN data available, you can generate sample data based on DMSP data:

```bash
julia conjunction_finder.jl --generate-sample
```

### Full Options

```bash
julia conjunction_finder.jl --dmsp path/to/dmsp_data.gz \
                          --elfin path/to/elfin_data.csv \
                          --mlt-threshold 0.5 \
                          --mlat-threshold 2.0 \
                          --min-duration 60 \
                          --output conjunction_events.csv
```

## Input Data Format

### DMSP Data
The DMSP data is expected in the binary SSJ format as described in the `read_ssj_file.pro` IDL code. The script can read both compressed (.gz) and uncompressed files.

### ELFIN Data
The ELFIN data should be in CSV format with the following columns:
- `timestamp`: Time in ISO format (e.g., "2023-01-01T00:00:00") or Unix timestamp
- `mlt`: Magnetic Local Time in hours (0-24)
- `mlat`: Magnetic Latitude in degrees (-90 to 90)

## Output

The script produces:
1. A CSV file with conjunction events (start time, end time, duration)
2. A JSON file with summary statistics
3. Plot images for each conjunction event showing MLT and MLAT values

## Example

```bash
# Download DMSP data and find conjunctions with generated ELFIN data
julia conjunction_finder.jl --dmsp https://www.ncei.noaa.gov/data/dmsp-space-weather-sensors/access/f17/ssj/2022/10/j5f1722288.gz --generate-sample --mlt-threshold 0.3 --mlat-threshold 1.5 --min-duration 120
```

This will:
1. Download DMSP data from the specified URL
2. Generate sample ELFIN data based on the DMSP data
3. Find conjunction events where MLT difference < 0.3 hours and MLAT difference < 1.5 degrees
4. Only include events lasting at least 120 seconds
5. Save results to conjunction_events.csv and conjunction_stats.json
6. Create plots for each conjunction event in the plots/ directory

## Project Structure

- `conjunction_finder.jl`: Main script for finding conjunctions
- `DMSPReader.jl`: Module for reading DMSP SSJ data files
- `ELFINReader.jl`: Module for reading and generating ELFIN data
- `run.sh`: Convenience script to run the conjunction finder
