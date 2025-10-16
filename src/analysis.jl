struct MLT{T} <: Real
    val::T
end

Base.promote_rule(::Type{MLT{T}}, ::Type{<:Real}) where {T} = MLT{T}
Base.:(+)(x::MLT, y::MLT) = MLT(mod(x.val + y.val, 24))

function parameters(df)
    return @transform(
        df,
        :log_A1 = (log10 ∘ A ∘ m1).(:model),
        :E_c1 = (E_c ∘ m1).(:model),
        :γ = (γ_len ∘ m1).(:model),
        :log_A2 = (log10 ∘ A ∘ m2).(:model),
        :E_c2 = (E_c ∘ m2).(:model),
        :κ = (κ ∘ m2).(:model),
        :Emin = Emin.(:model),
        :Δt = Δtrange.(:trange_elx, :trange_dmsp),
        :mlat, :maxAE, :mlt_elx, :mlt_dmsp
    )
end

function integrate_flux!(df, Emin)
    # Emin [keV]
    return @transform! df @astable begin
        :J1 = ((n_flux(Emin, 100) ∘ m1).(:model)) ./ 1000 # [1/cm^2/s/sr]
        :J2 = ((n_flux(Emin, 1000) ∘ m2).(:model)) ./ 1000
        :JE1 = ((e_flux(Emin, 100) ∘ m1).(:model)) ./ 1000 # [keV /cm^2/s/sr]
        :JE2 = ((e_flux(Emin, 1000) ∘ m2).(:model)) ./ 1000
        :J_e30 = ((n_flux(30, 1000) ∘ m2).(:model) .+ (n_flux(30, 1000) ∘ m1).(:model)) ./ 1000
        :JE_e30 = ((e_flux(30, 1000) ∘ m2).(:model) .+ (e_flux(30, 1000) ∘ m1).(:model)) ./ 1000
        :J_tot = :J1 + :J2
        :JE_tot = :JE1 + :JE2
        :R_J2 = :J2 ./ :J_tot
        :R_JE2 = :JE2 ./ :JE_tot
        :R_J_e30 = :J_e30 ./ :J_tot
        :R_JE_e30 = :JE_e30 ./ :JE_tot
        :E = :JE_tot ./ :J_tot
    end
end
