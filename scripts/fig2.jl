include(joinpath(@__DIR__, "data.jl"))

f = let df = combine(gdf, nrow => :n, renamecols = false) |> dropmissing!
    plt = data(df) * mapping(:mlt_bin => mlt2x => "MLT", :mlat_bin => "MLAT", :n) * visual(Heatmap) * mapping(col = :maxAE_bin)
    fg = draw(plt; figure = (; size = (1000, 300)), axis = m2axis)
    fg.grid[1].axis.yticks[] = 53:3:80
    easy_save("n_mlt_mlat")
end
