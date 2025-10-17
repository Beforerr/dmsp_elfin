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
            1e9, 5e8,         # J
            1e10, 8e9,        # JE
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

    @testset "Model Creation" begin
        # Create mock data
        df = DataFrame(
            mlat_bin = [65.0, 65.0, 70.0],
            mlt_bin = [6, 12, 6],
            ae_bin = ["[0, 100)", "[0, 100)", "[100, 300)"],
            n_samples = [50, 60, 45],
        )

        # Add parameter columns with statistics
        params = [:κ, :log_A1, :E_c1, :γ, :log_A2, :E_c2, :Emin, :J1, :J2, :JE1, :JE2]
        for param in params
            df[!, Symbol(string(param) * "_median")] = rand(3) .* 10
            df[!, Symbol(string(param) * "_q25")] = rand(3) .* 5
            df[!, Symbol(string(param) * "_q75")] = rand(3) .* 15
        end

        model = create_model_from_stats(df;
            mlat_bins=[65.0, 70.0],
            mlt_bins=[6, 12],
            ae_bins=["[0, 100)", "[100, 300)"],
            description="Test model",
            version="0.1.0"
        )

        @test model isa EmpiricalFluxModel
        @test model.version == "0.1.0"
        @test length(model.spatial_bins.mlat) == 2
        @test nrow(model.parameters) == 3
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

        model = create_model_from_stats(df;
            mlat_bins=[65.0, 70.0],
            mlt_bins=[6, 18],
            ae_bins=["[0, 100)"]
        )

        # Test nearest neighbor interpolation
        result = flux_parameters(model; mlat=65.5, mlt=7.0, ae=50.0, stat=:median)
        @test result isa FluxParameters
        @test result.n_samples == 50

        # Test different statistics
        result_q25 = flux_parameters(model; mlat=65.0, mlt=6.0, ae=50.0, stat=:q25)
        result_q75 = flux_parameters(model; mlat=65.0, mlt=6.0, ae=50.0, stat=:q75)

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

        model = create_model_from_stats(df;
            mlat_bins=[65.0],
            mlt_bins=[6],
            ae_bins=["[0, 100)"],
            description="Test display",
            version="1.0.0"
        )

        io = IOBuffer()
        show(io, model)
        output = String(take!(io))

        @test occursin("EmpiricalFluxModel", output)
        @test occursin("Test display", output)
        @test occursin("1.0.0", output)
    end

end

@testset "Integration Tests" begin
    @testset "End-to-End Model Usage" begin
        # Create a simple model
        df = DataFrame(
            mlat_bin = [65.0],
            mlt_bin = [6],
            ae_bin = ["[0, 100)"],
            n_samples = [100],
            κ_median = [3.5], κ_q25 = [3.0], κ_q75 = [4.0],
            log_A1_median = [10.0], log_A1_q25 = [9.5], log_A1_q75 = [10.5],
            E_c1_median = [5.0], E_c1_q25 = [4.5], E_c1_q75 = [5.5],
            γ_median = [-2.0], γ_q25 = [-2.5], γ_q75 = [-1.5],
            log_A2_median = [8.0], log_A2_q25 = [7.5], log_A2_q75 = [8.5],
            E_c2_median = [50.0], E_c2_q25 = [45.0], E_c2_q75 = [55.0],
            Emin_median = [30.0], Emin_q25 = [25.0], Emin_q75 = [35.0],
            J1_median = [1e9], J1_q25 = [8e8], J1_q75 = [1.2e9],
            J2_median = [5e8], J2_q25 = [4e8], J2_q75 = [6e8],
            JE1_median = [1e10], JE1_q25 = [8e9], JE1_q75 = [1.2e10],
            JE2_median = [8e9], JE2_q25 = [6e9], JE2_q75 = [1e10]
        )

        model = create_model_from_stats(df;
            mlat_bins=[65.0],
            mlt_bins=[6],
            ae_bins=["[0, 100)"]
        )

        # Test flux evaluation
        flux = evaluate_flux(model, 10.0; mlat=65.0, mlt=6.0, ae=50.0)
        @test flux > 0
        @test isfinite(flux)

        # Test spectrum
        energies = [1.0, 10.0, 100.0]
        spectrum = flux_spectrum(model, energies; mlat=65.0, mlt=6.0, ae=50.0)
        @test length(spectrum) == 3
        @test all(isfinite, spectrum)
        @test all(>(0), spectrum)

        # Test integrated fluxes
        J = number_flux(model, 1.0, 100.0; mlat=65.0, mlt=6.0, ae=50.0)
        JE = energy_flux(model, 1.0, 100.0; mlat=65.0, mlt=6.0, ae=50.0)
        @test J > 0 && isfinite(J)
        @test JE > 0 && isfinite(JE)

        # Test components
        comp = flux_components(model, 50.0; mlat=65.0, mlt=6.0, ae=50.0)
        @test comp.low_energy >= 0
        @test comp.high_energy >= 0
        @test comp.total ≈ comp.low_energy + comp.high_energy
    end
end
