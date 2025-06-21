import NaNMath as nm
export init_guess

function log_plec_model(E, p)
    A, γ, Ec = p
    log10_A = nm.log10(A)
    return log10_A .- γ .* nm.log10.(E) .- (E ./ Ec) .* nm.log10(ℯ)
end

function log_sbpl_model(E, p)
    A, γ1, γ2, Eb = p
    m = 1
    x = E ./ Eb
    return nm.log10(A) .- γ1 .* nm.log10.(x) .+ ((γ1 - γ2) ./ m) .* nm.log10.(1 .+ x .^ m)
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
    γ0 = 2.5
    Ec0 = E0 = energies[i]
    A = flux[i] / E0^γ0 * exp(E0 / Ec0)
    return [A, γ0, Ec0]
end

function init_guess(::typeof(log_sbpl_model), energies, flux)
    i = argmax(flux)
    m = 1
    γ1, γ2 = -5., 3.
    E0 = energies[i]
    Eb = 1e3
    A = flux[i] / exp(log_sbpl_model(E0, [1, γ1, γ2, Eb]))
    return [A, γ1, γ2, Eb]
end

function init_guess(::typeof(log_two_pop_model), energies, flux)
    i = argmax(flux)
    E0 = energies[i]
    Ec1 = E0
    Ec2 = 100E0
    γ1 = γ2 = 2.5
    A1 = flux[i] / E0^γ1 * exp(E0 / Ec1)
    A2 = A1 / 100
    return [A1, γ1, Ec1, A2, γ2, Ec2]
end

function init_guess(::typeof(log_plec_sbpl_model), energies, flux)
    p1 = [1e8, 3., 6e10]
    p2 = [8e7, -5.415023026784938, 2.9444485957997975, 63.69365352510846]
    return [p1..., p2...]
end

for model in (:log_plec_model, :log_sbpl_model, :log_two_pop_model, :log_plec_sbpl_model)
    sciml_model = Symbol(model, :_sciml)
    @eval $sciml_model(p, E) = $model(E, p)
    @eval export $model
    @eval export $sciml_model
    @eval init_guess(::typeof($sciml_model), energies, flux) = init_guess($model, energies, flux)
end