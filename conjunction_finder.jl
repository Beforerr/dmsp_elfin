#!/usr/bin/env julia

"""
DMSP-ELFIN Conjunction Finder

This script identifies conjunction events between DMSP and ELFIN satellites
based on their MLT (Magnetic Local Time) and MLAT (Magnetic Latitude) values.

A conjunction is defined when:
- |MLT_DMSP - MLT_ELFIN| < d1
- |MLAT_DMSP - MLAT_ELFIN| < d2
- This condition persists for at least T_min
"""

using ArgParse
using CSV
using DataFrames
using Dates
using Downloads
using GZip
using Interpolations
using JSON
using Plots
using Printf
using Statistics

using .ELFINReader



"""
    handle_mlt_wraparound(mlt_diff)

Handle the wraparound case for MLT differences (e.g., 23.5 - 0.5 = 23 should be 1).
"""
function handle_mlt_wraparound(mlt_diff)
    # MLT is in range [0, 24), so the maximum difference should be 12
    return min.(abs.(mlt_diff), 24 .- abs.(mlt_diff))
end

"""
    interpolate_data(dmsp_df::DataFrame, elfin_df::DataFrame, freq_seconds=1)

Interpolate data to a common time grid for comparison.
"""
function interpolate_data(dmsp_df::DataFrame, elfin_df::DataFrame, freq_seconds=1)
    # Find overlapping time range
    start_time = max(minimum(dmsp_df.timestamp), minimum(elfin_df.timestamp))
    end_time = min(maximum(dmsp_df.timestamp), maximum(elfin_df.timestamp))

    println("Interpolating data from $start_time to $end_time with frequency $(freq_seconds)s")

    # Create common time index
    common_times = start_time:Second(freq_seconds):end_time

    # Convert timestamps to float seconds for interpolation
    dmsp_times_sec = Float64.(Dates.value.(dmsp_df.timestamp .- DateTime(1970, 1, 1)) ./ 1000)
    elfin_times_sec = Float64.(Dates.value.(elfin_df.timestamp .- DateTime(1970, 1, 1)) ./ 1000)
    common_times_sec = Float64.(Dates.value.(common_times .- DateTime(1970, 1, 1)) ./ 1000)

    # Sort data by time to ensure proper interpolation
    dmsp_sorted_idx = sortperm(dmsp_times_sec)
    elfin_sorted_idx = sortperm(elfin_times_sec)

    # Create interpolation objects for DMSP data with sorted times
    dmsp_mlat_interp = linear_interpolation(
        dmsp_times_sec[dmsp_sorted_idx],
        dmsp_df.mlat[dmsp_sorted_idx],
        extrapolation_bc=Line()
    )
    dmsp_mlt_interp = linear_interpolation(
        dmsp_times_sec[dmsp_sorted_idx],
        dmsp_df.mlt[dmsp_sorted_idx],
        extrapolation_bc=Line()
    )

    # Create interpolation objects for ELFIN data with sorted times
    elfin_mlat_interp = linear_interpolation(
        elfin_times_sec[elfin_sorted_idx],
        elfin_df.mlat[elfin_sorted_idx],
        extrapolation_bc=Line()
    )
    elfin_mlt_interp = linear_interpolation(
        elfin_times_sec[elfin_sorted_idx],
        elfin_df.mlt[elfin_sorted_idx],
        extrapolation_bc=Line()
    )

    # Interpolate to common time grid
    dmsp_interp = DataFrame(
        timestamp=common_times,
        mlat=dmsp_mlat_interp.(common_times_sec),
        mlt=dmsp_mlt_interp.(common_times_sec)
    )

    elfin_interp = DataFrame(
        timestamp=common_times,
        mlat=elfin_mlat_interp.(common_times_sec),
        mlt=elfin_mlt_interp.(common_times_sec)
    )

    return dmsp_interp, elfin_interp
end

"""
    find_conjunctions(dmsp_df, elfin_df, mlt_threshold, mlat_threshold, min_duration_seconds)

Find conjunction events between DMSP and ELFIN satellites.
"""
function find_conjunctions(dmsp_df, elfin_df,
    mlt_threshold, mlat_threshold,
    min_duration_seconds)
    # Interpolate data to common time grid
    dmsp_interp, elfin_interp = interpolate_data(dmsp_df, elfin_df)

    # Calculate differences in MLT and MLAT
    mlt_diff = handle_mlt_wraparound(dmsp_interp.mlt .- elfin_interp.mlt)
    mlat_diff = abs.(dmsp_interp.mlat .- elfin_interp.mlat)

    # Create a mask for conjunction conditions
    conjunction_mask = (mlt_diff .< mlt_threshold) .& (mlat_diff .< mlat_threshold)

    # Find start and end indices of conjunction events
    state_changes = diff(vcat(0, conjunction_mask, 0))
    starts = findall(state_changes .== 1)
    ends = findall(state_changes .== -1) .- 1

    # Calculate duration of each event
    valid_events = []

    for i in eachindex(starts)
        start_idx = starts[i]
        end_idx = ends[i]
        start_time = dmsp_interp.timestamp[start_idx]
        end_time = dmsp_interp.timestamp[end_idx]
        duration = Dates.value(end_time - start_time) / 1000  # Duration in seconds

        if duration >= min_duration_seconds
            push!(valid_events, (
                start_time=start_time,
                end_time=end_time,
                duration_seconds=duration,
                start_idx=start_idx,
                end_idx=end_idx
            ))
        end
    end

    # Create DataFrame with conjunction events
    if !isempty(valid_events)
        events_df = DataFrame(valid_events)
        println("Found $(nrow(events_df)) conjunction events")
        return events_df, dmsp_interp, elfin_interp
    else
        println("No conjunction events found")
        return DataFrame(), dmsp_interp, elfin_interp
    end
end

"""
    analyze_conjunctions(events_df, dmsp_interp, elfin_interp)

Analyze conjunction events and generate summary statistics.
"""
function analyze_conjunctions(events_df, dmsp_interp, elfin_interp)
    if nrow(events_df) == 0
        return Dict("count" => 0)
    end

    # Calculate statistics
    durations = events_df.duration_seconds

    stats = Dict(
        "count" => nrow(events_df),
        "total_duration_seconds" => sum(durations),
        "min_duration_seconds" => minimum(durations),
        "max_duration_seconds" => maximum(durations),
        "mean_duration_seconds" => mean(durations),
        "median_duration_seconds" => median(durations)
    )

    # Plot each conjunction event
    for i in 1:nrow(events_df)
        plot_conjunction(events_df[i, :], dmsp_interp, elfin_interp)
    end

    return stats
end

"""
    save_results(events_df, stats, output_file="conjunction_events.csv", stats_file="conjunction_stats.json")

Save conjunction events and statistics to files.
"""
function save_results(events_df, stats, output_file="conjunction_events.csv", stats_file="conjunction_stats.json")
    if nrow(events_df) > 0
        CSV.write(output_file, events_df[:, [:start_time, :end_time, :duration_seconds]])
        println("Saved $(nrow(events_df)) conjunction events to $output_file")
    end

    # Save statistics as JSON
    open(stats_file, "w") do f
        print(f, JSON.json(stats, 2))
    end
    println("Saved statistics to $stats_file")

    # Print summary
    println("\nConjunction Analysis Summary:")
    println("Found $(stats["count"]) conjunction events")
    if stats["count"] > 0
        println("Total duration: $(round(stats["total_duration_seconds"]/60, digits=1)) minutes")
        println("Mean duration: $(round(stats["mean_duration_seconds"], digits=1)) seconds")
        println("Median duration: $(round(stats["median_duration_seconds"], digits=1)) seconds")
        println("Min duration: $(round(stats["min_duration_seconds"], digits=1)) seconds")
        println("Max duration: $(round(stats["max_duration_seconds"], digits=1)) seconds")
    end
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
    mlt_offset = 0.5 * sin.(range(0, 4π, length=n))
    elfin_df.mlt = mod.(elfin_df.mlt .+ mlt_offset, 24)

    # Add slight offset to MLAT (with some periods of close conjunction)
    mlat_offset = 2.0 * sin.(range(0, 6π, length=n))
    elfin_df.mlat = elfin_df.mlat .+ mlat_offset

    # Save to CSV
    CSV.write(output_file, elfin_df)
    println("Saved sample ELFIN data to $output_file")

    return elfin_df
end

"""
    parse_commandline()

Parse command line arguments.
"""
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--dmsp"
        help = "Path to DMSP data file or URL"
        default = "https://www.ncei.noaa.gov/data/dmsp-space-weather-sensors/access/f17/ssj/2022/10/j5f1722288.gz"
        "--elfin"
        help = "Path to ELFIN data file"
        "--mlt-threshold"
        help = "Maximum allowed difference in MLT (hours)"
        arg_type = Float64
        default = 0.5
        "--mlat-threshold"
        help = "Maximum allowed difference in MLAT (degrees)"
        arg_type = Float64
        default = 2.0
        "--min-duration"
        help = "Minimum duration for a conjunction event (seconds)"
        arg_type = Int
        default = 60
        "--output"
        help = "Output file for conjunction events"
        default = "conjunction_events.csv"
        "--generate-sample"
        help = "Generate sample ELFIN data for testing"
        action = :store_true
    end

    return parse_args(s)
end

"""
    main()

Main function to run the conjunction finder.
"""
function main()
    # Parse command line arguments
    args = parse_commandline()

    # Create data directory if it doesn't exist
    mkpath("data")

    # Handle DMSP data
    dmsp_file = args["dmsp"]
    if startswith(dmsp_file, "http")
        local_dmsp_file = joinpath("data", basename(dmsp_file))
        dmsp_file = DMSPReader.download_dmsp_data(dmsp_file, local_dmsp_file)
    end

    # Read DMSP data
    println("Reading DMSP data from $dmsp_file")
    dmsp_df = DMSPReader.read_dmsp_ssj_file(dmsp_file)
    println("Read $(nrow(dmsp_df)) DMSP records")

    # Handle ELFIN data
    elfin_file = args["elfin"]
    if args["generate_sample"]
        elfin_file = joinpath("data", "elfin_sample.csv")
        elfin_df = ELFINReader.generate_sample_elfin_data(dmsp_df, elfin_file)
    elseif isnothing(elfin_file)
        error("Please provide path to ELFIN data file or use --generate-sample")
    else
        println("Reading ELFIN data from $elfin_file")
        elfin_df = ELFINReader.read_elfin_data(elfin_file)
    end

    # Find conjunctions
    events_df, dmsp_interp, elfin_interp = find_conjunctions(
        dmsp_df, elfin_df,
        args["mlt_threshold"],
        args["mlat_threshold"],
        args["min_duration"]
    )

    # Analyze conjunctions
    stats = analyze_conjunctions(events_df, dmsp_interp, elfin_interp)

    # Save results
    save_results(events_df, stats, args["output"])
end