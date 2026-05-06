# DMSP-ELFIN Project

## Project layout

- `overleaf/` — manuscript (`index.tex`, `response.md`)
- `scripts/setup.jl` — common imports + `include` of src files
- `scripts/data.jl` — loads stats pipeline (`df` → `sdf` → `gdf`) into session scope; data cached as JLD2 via `produce_or_load`, loads fast from cache
- `scripts/fig1.jl` — ELFIN+DMSP conjunction demo plot; includes `setup.jl` directly (event-specific, no stats needed)
- `scripts/fig2.jl`–`fig5.jl` — statistical figure scripts; each includes `data.jl`