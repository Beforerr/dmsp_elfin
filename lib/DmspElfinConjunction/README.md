# DMSP-ELFIN Conjunction Finder [WIP]

This tool identifies conjunction events between DMSP and ELFIN satellites based on their MLT (Magnetic Local Time) and MLAT (Magnetic Latitude) values.

## What is a Conjunction Event?

A conjunction event is defined when:

- The absolute difference in MLT between DMSP and ELFIN is less than a threshold d₁
- The absolute difference in MLAT between DMSP and ELFIN is less than a threshold d₂
- This condition persists for at least Tₘᵢₙ seconds

