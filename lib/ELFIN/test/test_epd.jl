using Test
using ELFIN
using ELFIN: CDF


"""
Test to verify that our Julia EPD processing matches Python's implementation.
This test compares the precipitating flux calculation between Julia and Python.
"""

@testset "EPD Processing Verification" begin
    # Test parameters
    trange = ("2021-08-08", "2021-08-09")
    probe = "b"
    using Base.Iterators: product

    array_isapprox(x, y; kwargs...) = all(eachindex(x, y)) do i
        isnan(x[i]) && isnan(y[i]) && return true
        isapprox(x[i], y[i]; kwargs...)
    end

    no_update = true
    foreach(product((false, true), (:nflux, :eflux))) do (fullspin, type)
        py_results = @time ELFIN.py_epd(trange, probe; no_update, fullspin, type)
        jl_results = @time ELFIN.epd(trange, probe; no_update, fullspin, type)
        for (x, y) in zip(values(jl_results), values(py_results))
            @test array_isapprox(x, y; rtol = 1.0e-3)
        end
    end
end


@testset "Time" begin
    # Test parameters
    trange = ("2021-08-08", "2021-08-09")
    probe = "b"
    no_update = true
    using Base.Iterators: product
    cdf = ELFIN.load(trange, probe; no_update)
    var = "elb_pef_hs_time"
    @test CDF.tt2000_to_datetime_py(cdf.py[var]) == cdf[var]

    var = "elb_pef_hs_Epat_nflux"
end
