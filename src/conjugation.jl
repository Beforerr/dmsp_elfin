"""
    find_separate_conditions(continuous_ranges, elx_gei; kw...)

Find time ranges where Δmlt < Δmlt_max and Δmlat < Δmlat_max are satisfied separately.
"""
function find_separate_conditions(continuous_ranges, elx_gei; kw...)
    return filter(!isnothing, check_range_separate_conditions.(continuous_ranges, Ref(elx_gei); kw...))
end

"""
    find_simultaneous_conditions(continuous_ranges, elx_gei; kw...)

Find time ranges where both Δmlt < Δmlt_max and Δmlat < Δmlat_max are satisfied simultaneously.
"""
function find_simultaneous_conditions(continuous_ranges, elx_gei; kw...)
    return filter(!isnothing, check_range_conditions.(continuous_ranges, Ref(elx_gei); kw...))
end


function has_simultaneous_match(elx_mlt, elx_mlat, dmsp_mlt, dmsp_mlat; Δmlt_max = 1, Δmlat_max = 8)
    @inbounds for i in eachindex(elx_mlt, elx_mlat)
        e_mlt = elx_mlt[i]
        e_mlat = elx_mlat[i]
        for j in eachindex(dmsp_mlt, dmsp_mlat)
            if mlt_dist(e_mlt, dmsp_mlt[j]) < Δmlt_max && abs(e_mlat - dmsp_mlat[j]) < Δmlat_max
                return true
            end
        end
    end
    return false
end

"""
    check_range_conditions(ranges, elx_gei; ids = 16:18, Δt = Minute(10), kw...)

Check if a time range satisfies simultaneous MLT and MLAT conditions for DMSP satellites `ids`.

# Returns
- `(range=ranges, ids=valid_ids)` if conditions are met
- `nothing` if no DMSP satellite satisfies both conditions simultaneously
"""
function check_range_conditions(ranges, elx_gei; ids = 16:18, Δt = Minute(10), kw...)
    elx_mlt, elx_mlat = gei2mlt_mlat(tview(elx_gei, ranges))
    valid_ids = Int[]
    foreach(ids) do id
        dmsp_mlt, dmsp_mlat = get_mlt_mlat(extend(ranges, Δt), id)
        if has_simultaneous_match(elx_mlt, elx_mlat, dmsp_mlt, dmsp_mlat; kw...)
            push!(valid_ids, id)
        end
    end
    return !isempty(valid_ids) ? (range = ranges, ids = valid_ids) : nothing
end

"""
    check_range_separate_conditions(ranges, elx_gei; ids = 16:18, Δt = Minute(10), Δmlt_max = 1, Δmlat_max = 8)

Check if a time range satisfies MLT and MLAT conditions separately for DMSP satellites.
"""
function check_range_separate_conditions(ranges, elx_gei; ids = 16:18, Δt = Minute(10), Δmlt_max = 1, Δmlat_max = 8)
    elx_mlt, elx_mlat = gei2mlt_mlat(tview(elx_gei, ranges))
    valid_ids = Int[]
    foreach(ids) do id
        dmsp_mlt, dmsp_mlat = get_mlt_mlat(extend(ranges, Δt), id)
        if dist(mlt_dist, elx_mlt, dmsp_mlt) < Δmlt_max && dist(elx_mlat, dmsp_mlat) < Δmlat_max
            push!(valid_ids, id)
        end
    end
    return !isempty(valid_ids) ? (range = ranges, ids = valid_ids) : nothing
end
