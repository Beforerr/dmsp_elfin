include(joinpath(@__DIR__, "setup.jl"))

using SpectralModels: A, κ, E_c, n_flux, e_flux

const γ_len = @optic _.γ
const m1 = @optic _.model1
const m2 = @optic _.model2
const Emin = @optic _.Emin

tr = Date("2020-01-01"), Date("2022-10-01")
ids = 16:18
df_a = produce(tsplit(tr, Month), "a", ids)
df_b = produce(tsplit(tr, Month), "b", ids)

df_a.elfin .= :a
df_b.elfin .= :b
df = vcat(df_a, df_b)

# Results without sanitizing DMSP flux for comparison
# df2 = @rtransform! df $AsTable = fit_two_flux(:flux_dmsp, :flux_elx; threshold = FLUX_THRESHOLD)

sdf = let Emin = 0.03
    @chain parameters(df) begin
        @rsubset!(20 > :κ > 1, :E_c1 < 20, :E_c2 < 20, abs(:γ) < 20, 0 < :log_A1 < 30, 0 < :log_A2 < 20)
        integrate_flux!(Emin)
    end
end

bin(x, δx) = (x + δx / 2) - (x + δx / 2) % δx
demlt(x::T) where {T} = x == 24 ? zero(T) : x
fmt(from, to, i; leftclosed, rightclosed) = from + (to - from) / 2

cut_AE(df; breaks = [0, 100, 300]) = @transform(df, :maxAE_bin = cut(:maxAE, breaks; extend = true))

make_gdf(df; mlat_binedges = 52.5:1:80.5, Δmlt = 1) = @chain df begin
    cut_AE
    @transform(
        :mlt_bin = demlt.(bin.(mlt_mean.(:mlt_elx, :mlt_dmsp), Δmlt)),
        :mlat_bin = Array(cut(abs.(:mlat), mlat_binedges; extend = missing, labels = fmt)),
    )
    @groupby(:mlt_bin, :mlat_bin, :maxAE_bin)
end

gdf = make_gdf(sdf)

_make_tdf(gdf, func, vars; n = 6) = @chain gdf begin
    combine(_variable.(vars) .=> func, nrow => :n, renamecols = false)
    dropmissing!
    @rsubset(:n >= n)
end

const m2axis = (; xtickformat = mlt_tickformat, xticks = 0:3:24)
