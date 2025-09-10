function gei(probe, t0, t1)
    probe = "el$(probe)"
    url_pattern = "$BASE_URL/$probe/l1/state/defn/{Y}/$(probe)_l1_state_defn_{Y}{M:02d}{D:02d}_v\\d+.cdf"
    var = get_product(
        url_pattern, "$(probe)_pos_gei", t0, t1;
        use_file_list = true,
    )
    return DimArray(var)
end

function gei(trange, probe)
    return gei(probe, trange[1], trange[2])
end
