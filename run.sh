#!/bin/bash

# DMSP-ELFIN Conjunction Finder Runner
# This script downloads DMSP data and runs the conjunction finder

# Create data directory if it doesn't exist
mkdir -p data

# Default values
DMSP_URL="https://www.ncei.noaa.gov/data/dmsp-space-weather-sensors/access/f17/ssj/2022/10/j5f1722288.gz"
MLT_THRESHOLD=0.5
MLAT_THRESHOLD=2.0
MIN_DURATION=60

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --dmsp)
      DMSP_URL="$2"
      shift 2
      ;;
    --mlt-threshold)
      MLT_THRESHOLD="$2"
      shift 2
      ;;
    --mlat-threshold)
      MLAT_THRESHOLD="$2"
      shift 2
      ;;
    --min-duration)
      MIN_DURATION="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

echo "Running DMSP-ELFIN Conjunction Finder"
echo "DMSP URL: $DMSP_URL"
echo "MLT Threshold: $MLT_THRESHOLD"
echo "MLAT Threshold: $MLAT_THRESHOLD"
echo "Min Duration: $MIN_DURATION seconds"

# Run the conjunction finder
julia conjunction_finder.jl \
  --dmsp "$DMSP_URL" \
  --generate-sample \
  --mlt-threshold "$MLT_THRESHOLD" \
  --mlat-threshold "$MLAT_THRESHOLD" \
  --min-duration "$MIN_DURATION"

echo "Conjunction finder completed"
