include(joinpath(@__DIR__, "data.jl"))

# Transform dataframe with range categories
let params = [:R_JE_e30 => L"$J_E^{>30 \text{keV}} / J_E$", :E => L"$\bar{E}$ [keV]", :κ => L"$κ$"]
    labels = last.(params)
    pdf_df = @chain cut_AE(sdf) begin
        @rtransform(
            :MLT = ifelse(12 <= :mlt_elx <= 21, "12-21", :mlt_elx < 9 ? "0-9" : missing),
            :MLAT = 60 <= abs(:mlat) < 65 ? "60-65" : (abs(:mlat) >= 65 ? ">65" : missing)
        )
        dropmissing([:maxAE_bin, :MLT, :MLAT])
    end

    plt = data(pdf_df) *
        mapping(params, color = :MLAT, linestyle = :MLT) *
        mapping(col = :maxAE_bin, row = AlgebraOfGraphics.dims(1) => renamer(labels)) *
        visual(Lines; linewidth = 2) *
        AlgebraOfGraphics.density(; datalimits = x -> quantile(x, [0.04, 0.96]))

    fg = draw(plt; facet = (; linkxaxes = :none, linkyaxes = :minimal), axis = (xlabel = "", yscale = log10))

    Makie.xlims!.(fg.grid[1, :], -0.05, 0.9)
    Makie.xlims!.(fg.grid[2, :], -0.2, 12)
    Makie.ylims!.(fg.grid[2, :], 10^(-2.4), 10^(-0.2))
    Makie.xlims!.(fg.grid[3, :], 1, 7.8)
    for (i, j) in Iterators.product(1:3, 1:3)
        lbl = "($(('a':'z')[i])$j)"
        text!(fg.grid[i, j].axis, 0.92, 0.97; text = lbl, space = :relative, align = (:right, :top), font = :bold)
    end
    easy_save("pdf_params_ae_mlt_mlat")
end
