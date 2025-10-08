using DmspElfinConjunction: bin_mlat_times

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


"""
    find_matched_mlat_conditions(continuous_ranges, elx_gei; kw...)

Find time ranges where corresponding MLT differences are also within threshold after matching binned MLAT.
"""
function find_matched_mlat_conditions(continuous_ranges, elx_gei; ids = 16:18, Δt = Minute(10), kw...)
    return filter(!isnothing, check_matched_mlat_Δmlt.(continuous_ranges, Ref(elx_gei); ids, Δt, kw...))
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
        dmsp_mlt, dmsp_mlat = get_mlt_mlat(id, extend(ranges, Δt))
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
        dmsp_mlt, dmsp_mlat = get_mlt_mlat(id, extend(ranges, Δt))
        if dist(mlt_dist, elx_mlt, dmsp_mlt) < Δmlt_max && dist(elx_mlat, dmsp_mlat) < Δmlat_max
            push!(valid_ids, id)
        end
    end
    return !isempty(valid_ids) ? (range = ranges, ids = valid_ids) : nothing
end

"""
    check_matched_mlat_Δmlt(elx_mlt, elx_mlat, dmsp_mlt, dmsp_mlat; δmlat = 0.5, Δmlt_max = 1, δt = Millisecond(1000))

First bin MLAT data, then check if MLT differences within overlapping bins are within Δmlt_max.


# Returns
`true` if any overlapping MLAT bins have MLT difference < Δmlt_max, `false` otherwise.
"""
function check_matched_mlat_Δmlt(elx_mlt, elx_mlat, dmsp_mlt, dmsp_mlat; δmlat = 0.5, Δmlt_max = 1, δt = Millisecond(1000))
    elx_bins = bin_mlat_times(elx_mlat; δmlat, δt)
    dmsp_bins = bin_mlat_times(dmsp_mlat; δmlat, δt)
    common_mlat_bins = intersect(keys(elx_bins), keys(dmsp_bins))
    return any(common_mlat_bins) do mlat_bin
        # Get MLT values within these time ranges
        elx_tranges = elx_bins[mlat_bin]
        dmsp_tranges = dmsp_bins[mlat_bin]
        any(Base.Iterators.product(elx_tranges, dmsp_tranges)) do (elx_trange, dmsp_trange)
            mlt_dist(
                local_mlt_mean(tview(elx_mlt, elx_trange)), # todo what is the mean of [1,23]
                local_mlt_mean(tview(dmsp_mlt, dmsp_trange))
            ) < Δmlt_max
        end
    end
end

function check_matched_mlat_Δmlt(ranges, elx_gei; ids = 16:18, Δt = Minute(10), Δmlt_max = 1)
    _elx_gei = tview(elx_gei, ranges)
    isempty(_elx_gei) && return nothing
    elx_mlt, elx_mlat = gei2mlt_mlat(_elx_gei)
    all(isnan, elx_mlat) && return nothing
    valid_ids = Int[]
    foreach(ids) do id
        dmsp_mlt, dmsp_mlat = get_mlt_mlat(id, extend(ranges, Δt))
        try
            if check_matched_mlat_Δmlt(elx_mlt, elx_mlat, dmsp_mlt, dmsp_mlat; Δmlt_max)
                push!(valid_ids, id)
            end
        catch
            Main.@autoinfiltrate
            @warn "Failed to check matched MLAT conditions for $id with ranges $(ranges)"
            rethrow()
        end
    end
    return !isempty(valid_ids) ? (range = ranges, ids = valid_ids) : nothing
end
