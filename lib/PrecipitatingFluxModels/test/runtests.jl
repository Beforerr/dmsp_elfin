using Test
using PrecipitatingFluxModels
using DataFrames

@testset "PrecipitatingFluxModels.jl" begin

    @testset "AE Bin Selection" begin
        ae_bins = ["[0, 100)", "[100, 300)", "[300, Inf)"]

        @test get_ae_bin(50.0, ae_bins) == "[0, 100)"
        @test get_ae_bin(150.0, ae_bins) == "[100, 300)"
        @test get_ae_bin(400.0, ae_bins) == "[300, Inf)"
        @test get_ae_bin(0.0, ae_bins) == "[0, 100)"
        @test get_ae_bin(99.9, ae_bins) == "[0, 100)"
        @test get_ae_bin(100.0, ae_bins) == "[100, 300)"
    end

    @testset "FluxParameters" begin
        params = FluxParameters(
            10.0, 5.0, -2.0,  # ExpPow
            8.0, 50.0, 3.5,   # Kappa
            30.0,             # Emin
            1.0e9, 5.0e8,         # J
            1.0e10, 8.0e9,        # JE
            65.0, 6.0, 150.0, 100  # metadata
        )

        @test params.log_A1 == 10.0
        @test params.κ == 3.5
        @test params.mlat == 65.0

        # Test conversion to spectral model
        model = to_spectral_model(params)
        @test model isa TwoStepModel
        @test model.Emin == 30.0
    end

    @testset "Parameter Interpolation" begin
        # Create simple test model
        df = DataFrame(
            mlat_bin = [65.0, 65.0, 70.0, 70.0],
            mlt_bin = [6, 18, 6, 18],
            ae_bin = repeat(["[0, 100)"], 4),
            n_samples = [50, 50, 50, 50],
        )

        params = [:κ, :log_A1, :E_c1, :γ, :log_A2, :E_c2, :Emin, :J1, :J2, :JE1, :JE2]
        for param in params
            df[!, Symbol(string(param) * "_median")] = [1.0, 2.0, 3.0, 4.0]
            df[!, Symbol(string(param) * "_q25")] = [0.5, 1.5, 2.5, 3.5]
            df[!, Symbol(string(param) * "_q75")] = [1.5, 2.5, 3.5, 4.5]
        end

        model = nothing

        # Test nearest neighbor interpolation
        result = flux_parameters(model; mlat = 65.5, mlt = 7.0, ae = 50.0, stat = :median)
        @test result isa FluxParameters
        @test result.n_samples == 50

        # Test different statistics
        result_q25 = flux_parameters(model; mlat = 65.0, mlt = 6.0, ae = 50.0, stat = :q25)
        result_q75 = flux_parameters(model; mlat = 65.0, mlt = 6.0, ae = 50.0, stat = :q75)

        @test result_q25.κ < result_q75.κ
    end

    @testset "Model Type Display" begin
        df = DataFrame(
            mlat_bin = [65.0],
            mlt_bin = [6],
            ae_bin = ["[0, 100)"],
            n_samples = [50],
        )

        params = [:κ, :log_A1, :E_c1, :γ, :log_A2, :E_c2, :Emin, :J1, :J2, :JE1, :JE2]
        for param in params
            df[!, Symbol(string(param) * "_median")] = [1.0]
            df[!, Symbol(string(param) * "_q25")] = [0.5]
            df[!, Symbol(string(param) * "_q75")] = [1.5]
        end

        model = nothing

        io = IOBuffer()
        show(io, model)
        output = String(take!(io))

        @test occursin("EmpiricalFluxModel", output)
        @test occursin("Test display", output)
        @test occursin("1.0.0", output)
    end

end
