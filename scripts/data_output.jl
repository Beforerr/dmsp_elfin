include("../scripts/data.jl")

unique(df.mechanism), counts(Int.(df.mechanism))

tranges = @chain df begin
    @groupby(:elfin, :id)
    @combine $AsTable = (trs = find_continuous_timeranges(sort!(DateTime.(first.(:trange_elx))), Minute(2)); (t0 = first.(trs), t1 = last.(trs)))
    @orderby(:elfin, :id, :t0)
end

@info tranges_summary = @chain tranges begin
    @groupby(:elfin, :id)
    combine(nrow => :n, :t0 => minimum => :t0_min, :t1 => maximum => :t1_max)
    sort([:elfin, :id])
end

## Output the dataset as a CSV file
using CSV

let colnames = [:elfin, :id, :mlat, :mlt_elx, :mlt_dmsp, :maxAE, :trange_elx, :trange_dmsp, :flux_elx, :flux_dmsp, :model]
    select!(df, colnames, Not(colnames))
    CSV.write(joinpath(@__DIR__, "../files/events.csv"), df)
end
