function log_plec_model(E, p)
    A, γ, Ec = p
    return nm.log(A) .- γ .* nm.log.(E) .- (E ./ Ec)
end

function log_plec_sbpl_model(E, p)
    A1, γ1, Ec1,    # PLEC params
        A2, γ2, γ3, Eb = p  # SBPL params
    m = 1

    # Component 1: PLEC
    f1 = @. A1 * E^(-γ1) * exp(-E / Ec1)

    # Component 2: SBPL
    x = E ./ Eb
    f2 = @. A2 * x^(-γ2) * (1 + x^m)^((γ2 - γ3) / m)

    return log10.(f1 .+ f2)
end

function log_two_pop_model(p, E)
    A1, γ1, Ec1, A2, γ2, Ec2 = p
    f1 = A1 .* E .^ (-γ1) .* exp.(-E ./ Ec1)
    f2 = A2 .* E .^ (-γ2) .* exp.(-E ./ Ec2)
    flux = f1 .+ f2
    return NaNMath.log10.(flux)
end

function init_guess(::typeof(log_plec_model), energies, flux)
    i = argmax(flux)
    γ0 = 0.5
    E0 = energies[i]
    Ec0 = energies[end]
    A = flux[i] / E0^γ0 * exp(E0 / Ec0)
    return [A, γ0, Ec0]
end

function init_guess(::typeof(log_sbpl_model), energies, flux)
    i = argmax(flux)
    m = 1
    γ1, γ2 = -5.0, 3.0
    E0 = energies[i]
    Eb = 1.0e3
    A = flux[i] / exp(log_sbpl_model(E0, [1, γ1, γ2, Eb]))
    return [A, γ1, γ2, Eb]
end

function init_guess(::typeof(log_two_pop_model), energies, flux)
    i = argmax(flux)
    E0 = energies[i]
    Ec1 = E0
    Ec2 = 100.0e0
    γ1 = γ2 = 2.5
    A1 = flux[i] / E0^γ1 * exp(E0 / Ec1)
    A2 = A1 / 100
    return [A1, γ1, Ec1, A2, γ2, Ec2]
end