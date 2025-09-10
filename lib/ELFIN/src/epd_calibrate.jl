"""
    epd_calibrate(raw_counts; dead_time_correction=true)

Apply calibration to raw EPD counts to convert to physical flux units.

# Arguments
- `raw_counts`: Raw count data from EPD
- `dead_time_correction`: Apply dead time correction (default: true)

# Returns
- Calibrated flux data
"""
function epd_calibrate(raw_counts; dead_time_correction = true)
    # This is a simplified calibration - in practice, this would use
    # detailed calibration coefficients from instrument characterization

    calibrated = copy(raw_counts)

    if dead_time_correction
        # Apply dead time correction (simplified)
        dead_time = 1.5e-6  # 1.5 μs typical for EPD
        calibrated = calibrated ./ (1 .- calibrated .* dead_time)
    end

    # Apply energy-dependent geometric factors and efficiency corrections
    # (This would normally use detailed calibration files)
    geom_factors = fill(0.031256, length(EPD_ENERGY_BINS))  # Simplified geometric factor

    for (i, gf) in enumerate(geom_factors)
        if ndims(calibrated) == 3  # time × energy × sector
            calibrated[:, i, :] ./= gf
        elseif ndims(calibrated) == 2  # time × energy
            calibrated[:, i] ./= gf
        end
    end

    return calibrated
end
