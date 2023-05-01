"""Eich scaling(NF 53 093031) & B. Sieglin PPCF 55 (2013) 124039"""
function λq_sieglin(R, P_SOL, Bpol, νₑ)
    # Eich scaling(NF 53 093031)
    λ_q = 1.35 * (P_SOL)^(-0.02) * R^0.04 * Bpol^(-0.92) * νₑ^0.42 * 1.e-3
    #From B. Sieglin PPCF 55 (2013) 124039
    S = 0.09 * (νₑ / 1.e19)^1.02 * Bpol^(-1.01) * (R / 3.0) * 1.e-3

    return λ_q + 1.64 * S
end

"""Eich scaling (NF 53 093031)"""
λq_eich(R::T, P_SOL::T, Bpol::T, ϵ::T) where {T<:Real} = 1.35e-3 * P_SOL^(-0.02) * R^0.04 * Bpol^(-0.92) * ϵ^(0.42) # Eich scaling (NF 53 093031)

λq_eich(plasma::LengyelModelPlasmaParameters) = λq_eich(plasma.R_omp, plasma.P_SOL, plasma.Bpol_omp, plasma.ϵ)

λq_loarte(P_SOL::T, B0::T, q95iter::T) where {T<:Real} = 0.00265 * P_SOL^0.38 * B0^(-0.71) * q95iter^0.3