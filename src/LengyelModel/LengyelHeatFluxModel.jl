import FuseUtils
import SimulationParameters
import ADAS
import IMASDD
using Formatting

mutable struct LengyelModelResults{T}
    q_poloidal_omp::T
    q_parallel_omp::T
    q_rad::T
    q_parallel_target_unprojected::T
    q_parallel_target_projected::T
    q_perp_target::T
    q_perp_target_spread::T
    zeff_up::T
    λ_target::T
    f_pol_projection::T
    f_perp_projection::T
end

function Base.show(io::IO, r::LengyelModelResults)
    for f in propertynames(r)
        println(io, "$f: $(getfield(r,f))")
    end
end

LengyelModelResults() = LengyelModelResults((0.0 for d in fieldnames(LengyelModelResults))...)

mutable struct LengyelModel <: DivertorHeatFluxModel
    parameters::LengyelModelParameters
    results::LengyelModelResults
end

LengyelModel() = LengyelModel(LengyelModelParameters{Float64}())

function LengyelModel(par::LengyelModelParameters)
    SimulationParameters.setup_parameters!(par)
    LengyelModel(par, LengyelModelResults())
end

(model::LengyelModel)() = model.results = compute_lengyel_model(model.parameters)

""" Perform weighted cooling rate integral over specified temperature interval
# Inputs:  Tmin   minimum temperature for integral (eV)
#          Tmax   maximum temperature for integral (eV)
#          Texp   Exponent for temperature weighting
#          Lexp   Exponent for cooling rate weighting
#          Zimp   Z of impurity (Ar: 18; Kr: 36; Xe: 54)
"""
function V_legyel_ADAS(Tmin::Float64, Tmax::Float64, f_imp::Float64, imp::Union{String,Symbol}; N::Int64=500, ne::Float64=1e20, Zeff_exp::Float64=-0.3, Texp::Float64=0.5, Lexp::Float64=1.0, κ0=2390.0)
    data = ADAS.get_cooling_rates(imp)
    zeff = ADAS.get_Zeff(imp)
    Lz = data.Lztot
    T = collect(LinRange(Tmin, Tmax, N))
    int = [T_ .^ Texp .* zeff(f_imp, ne, T_) .^ (Zeff_exp) .* Lz(ne, T_) .^ Lexp for T_ in T]
    return sqrt.(FuseUtils.trapz(T, int) * f_imp * κ0 * 2)
end

function V_legyel_ADAS(Tmin::Float64, Tmax::Float64, f_imps::Vector{Float64}, imps::Vector{<:Union{String,Symbol}}; ne::Float64=1e20, Zeff_exp::Float64=-0.3, Texp::Float64=0.5, Lexp::Float64=1.0, κ0=2390.0, N::Int64=500)
    zeff = ADAS.get_Zeff(imps)
    T = collect(LinRange(Tmin, Tmax, N))
    int = 0.0
    for (f_imp, imp) in zip(f_imps, imps)
        data = ADAS.get_cooling_rates(imp)
        int += sqrt.(FuseUtils.trapz(T, [T_ .^ Texp .* zeff(f_imps, ne, T_) .^ (Zeff_exp) .* data.Lztot(ne, T_) .^ Lexp for T_ in T]) * f_imp * κ0 * 2.0)
    end
    return int
end

V_legyel_ADAS(s, i) = V_legyel_ADAS(i.T_down, s.T_up, s.f_imp, s.imp; ne=s.n_up, Zeff_exp=i.Zeff_exp, Texp=i.Texp, Lexp=i.Lexp, κ0=i.κ0)

function compute_lengyel_model(par::LengyelModelParameters)
    r = LengyelModelResults()
    r.λ_target = par.sol.λ_omp * par.target.f_omp2target_expansion
    r.f_pol_projection = tan(par.target.α_sp)
    r.f_perp_projection = 1.0 / sin(par.target.θ_sp)
    r.q_poloidal_omp = compute_q_poloidal_omp(par)
    r.q_parallel_omp = compute_q_parallel_omp(par)
    r.q_rad = compute_qrad(par)
    r.q_parallel_target_unprojected = sqrt(max(0.0, r.q_parallel_omp^2.0 - r.q_rad^2.0))
    r.q_parallel_target_projected = r.q_parallel_target_unprojected / par.target.f_omp2target_expansion * (par.plasma.R_omp  / par.target.R)
    r.q_perp_target = r.q_parallel_target_projected * r.f_pol_projection
    r.q_perp_target_spread = r.q_perp_target / par.target.f_spread_pfr
    r.zeff_up = compute_zeff_up(par)
    return r
end

function compute_heat_channel_area(R::T, λ_q::T) where {T<:Float64}
    @assert (R > 0.0 && λ_q > 0)
    return 2π * R * λ_q
end

compute_q_parallel_omp(p::LengyelModelParameters) = compute_q_parallel_omp(p.plasma.P_SOL, p.plasma.R_omp, p.sol.λ_omp, p.plasma.Bpol_omp, p.plasma.Bt_omp)

compute_q_poloidal_omp(p::LengyelModelParameters) = compute_q_poloidal_omp(p.plasma.P_SOL, p.plasma.R_omp, p.sol.λ_omp)

function compute_q_poloidal_omp(P_SOL::T, R::T, λ_q::T) where {T<:Float64}
    @assert (P_SOL >= 0.0)
    return P_SOL / compute_heat_channel_area(R, λ_q)
end

function compute_q_parallel_omp(P_SOL::T, R::T, λ_q::T, Bpol::T, Bt::T) where {T<:Float64}
    @assert (P_SOL >= 0.0 && Bpol > 0 && Bt > 0.0)
    return P_SOL / compute_heat_channel_area(R, λ_q) / sin(atan(Bpol / Bt))
end

compute_qrad(p::LengyelModelParameters) = compute_qrad(p.sol, p.integral)

compute_qrad(s::LengyelModelSOLParameters, i::LengyelIntegralParameters) = s.f_adhoc * s.n_up * s.T_up * V_legyel_ADAS(s, i)

function compute_zeff_up(par::LengyelModelParameters)
    zeff = ADAS.get_Zeff(par.sol.imp)
    return zeff(par.sol.f_imp, par.sol.n_up, par.sol.T_up)
end

"""
    show_summary(model::LengyelModel)

Print summary of LengyelModel setup and simulation results
"""
function show_summary(model::LengyelModel)
    p = model.parameters
    r = model.results
    printfmtln("Upstream")
    printfmtln("├─ {:<22} = {:.2f} MW", "P_SOL", p.plasma.P_SOL / 1e6)
    printfmtln("├─ {:<22} = {:.2f} eV", "Te_up", p.sol.T_up)
    printfmtln("├─ {:<22} = {:.2e} m⁻³", "ne_up", p.sol.n_up)
    printfmtln("├─ {:<22} = {:.1f} T", "Bp_omp", p.plasma.Bpol_omp)
    printfmtln("├─ {:<22} = {:.1f} T", "Bt_omp", p.plasma.Bt_omp)
    printfmtln("├─ {:<22} = {:.1f} m", "R_omp", p.plasma.R_omp)
    printfmtln("├─ {:<22} = {:.1f} MW.T/m", "Psol⋅Bp/R", p.plasma.P_SOL / 1e6 / p.plasma.R_omp * p.plasma.Bpol_omp)
    printfmtln("├─ {:<22} = {:.4f} m", "λ_omp", p.sol.λ_omp)
    printfmtln("└─ {:<22} = {:.4f} m²", "area_omp", compute_heat_channel_area(p.plasma.R_omp, p.sol.λ_omp))
    printfmtln("")
    printfmtln("Target")
    printfmtln("├─ {:<22} = {:.4f}", "f_omp2target_expansion", p.target.f_omp2target_expansion)
    printfmtln("├─ {:<22} = {:.4f}", "λ_target", r.λ_target)
    printfmtln("├─ {:<22} = {:.4f} m²", "area_target", compute_heat_channel_area(p.target.R, r.λ_target))
    printfmtln("├─ {:<22} = {:.1f} deg", "α_pitch", p.target.α_sp * 180 / π)
    printfmtln("├─ {:<22} = {:.1f}", "f_pol_projection", r.f_pol_projection)
    printfmtln("├─ {:<22} = {:.1f} deg", "θ_target", p.target.θ_sp * 180 / π)
    printfmtln("└─ {:<22} = {:.1f}", "f_perp_projection", r.f_perp_projection)
    printfmtln("")
    printfmtln("Transport")
    printfmtln("└─ {:<22} = {:.1f}", "spread_pfr", p.target.f_spread_pfr)
    printfmtln("")
    printfmtln("Impurities")
    for (i, f) in zip(p.sol.imp, p.sol.f_imp)
        printfmtln("├─ {:<22} = {:3.3f}%  ", string(i), f * 100)
    end
    printfmtln("└─ {:<22} = {:.2f} ", "Zeff_up", r.zeff_up)
    printfmtln("")
    printfmtln("Lengyel model output")
    printfmtln("├─ {:<22} = {:.2f} MW/m^2", "q_poloidal_omp", r.q_poloidal_omp / 1e6)
    printfmtln("├─ {:<22} = {:.2f} MW/m^2", "q_parallel_omp", r.q_parallel_omp / 1e6)
    printfmtln("├─ {:<22} = {:.2f} MW/m^2 ", "q_rad", (r.q_rad) / 1e6)
    printfmtln("├─ {:<22} = {:.2f} MW/m^2 ", "q_rad_effective", (r.q_parallel_omp - r.q_parallel_target_unprojected) / 1e6)
    printfmtln("├─ {:<22} = {:.2f} MW/m^2", "q_parallel_target", r.q_parallel_target_unprojected / 1e6)
    printfmtln("├─ {:<22} = {:.2f} MW/m^2 ", "q_perp_target", r.q_perp_target / 1e6)
    printfmtln("└─ {:<22} = {:.2f} MW/m^2 ", "q_perp_target_spread", r.q_perp_target_spread / 1e6)
end
