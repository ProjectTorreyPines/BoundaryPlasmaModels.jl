module StangebyHeatFluxModel
import SimulationParameters
import ADAS
using IMAS
using Format
using NLsolve
import ..DivertorHeatFluxModel
import ..DivertorHeatFluxModelParameters

export StangebyModelParameters, StangebyModel

include("Stangeby_parameters.jl")

mutable struct StangebyModelResults{T}
    q_poloidal_omp::T
    q_parallel_omp::T
    q_rad::T
    q_parallel_target_unprojected::T
    q_parallel_target_projected::T
    q_perp_target::T
    q_perp_target_spread::T
    λ_target::T
    f_pol_projection::T
    f_perp_projection::T
    T_target::T
end

function Base.show(io::IO, r::StangebyModelResults)
    for f in propertynames(r)
        println(io, "$f: $(getfield(r,f))")
    end
end

StangebyModelResults() = StangebyModelResults((0.0 for d in fieldnames(StangebyModelResults))...)

mutable struct StangebyModel <: DivertorHeatFluxModel
    parameters::StangebyModelParameters
    results::StangebyModelResults
end

StangebyModel() = StangebyModel(StangebyModelParameters{Float64}())

function StangebyModel(par::StangebyModelParameters)
    SimulationParameters.setup_parameters!(par)
    StangebyModel(par, StangebyModelResults())
end

(model::StangebyModel)() = model.results = compute_Stangeby_model(model.parameters)

function compute_Stangeby_model(par::StangebyModelParameters)
    r = StangebyModelResults()
    r.λ_target = par.sol.λ_omp * par.target.f_omp2target_expansion
    r.f_pol_projection = tan(par.target.α_sp)
    r.f_perp_projection = 1.0 / sin(par.target.θ_sp)
    r.q_poloidal_omp = compute_q_poloidal_omp(par)
    r.q_parallel_omp = compute_q_parallel_omp(par)
    r.q_rad = compute_qrad(par, r)
    r.q_parallel_target_unprojected = sqrt(max(0.0, r.q_parallel_omp^2.0 - r.q_rad^2.0))
    r.q_parallel_target_projected = r.q_parallel_target_unprojected / par.target.f_omp2target_expansion * (par.plasma.R_omp  / par.target.R)
    r.q_perp_target = r.q_parallel_target_projected * r.f_pol_projection
    r.q_perp_target_spread = r.q_perp_target / par.target.f_spread_pfr
    return r
end

function compute_heat_channel_area(R::T, λ_q::T) where {T<:Float64}
    @assert (R > 0.0 && λ_q > 0)
    return 2π * R * λ_q
end

compute_q_parallel_omp(p::StangebyModelParameters) = compute_q_parallel_omp(p.plasma.P_SOL, p.plasma.R_omp, p.sol.λ_omp, p.plasma.Bpol_omp, p.plasma.Bt_omp)

compute_q_poloidal_omp(p::StangebyModelParameters) = compute_q_poloidal_omp(p.plasma.P_SOL, p.plasma.R_omp, p.sol.λ_omp)

function compute_q_poloidal_omp(P_SOL::T, R::T, λ_q::T) where {T<:Float64}
    @assert (P_SOL >= 0.0)
    return P_SOL / compute_heat_channel_area(R, λ_q)
end

function compute_q_parallel_omp(P_SOL::T, R::T, λ_q::T, Bpol::T, Bt::T) where {T<:Float64}
    @assert (P_SOL >= 0.0 && Bpol > 0 && Bt > 0.0)
    return P_SOL / compute_heat_channel_area(R, λ_q) / sin(atan(Bpol / Bt))
end






function get_Te_target(Tᵤ, q_parallel_omp, L_para, f_cond, κ_0; Te_init=10.0)
    G(Tₑₜ) = @. abs(Tᵤ)^(7 / 2) - (abs(Tₑₜ[1])^(7 / 2) + 7 / 2 * q_parallel_omp * L_para * f_cond / κ_0 / (1.0 - fcool(Tₑₜ[1])))
    sol = nlsolve(G, [Te_init])
    return sol.zero[1]
end

# from https://doi.org/10.1088/1361-6587/ad2b90 (J. H. Nichols et al 2024 Plasma Phys. Control. Fusion 66 045013)
a = [-1.171, 2.101, -3.212, 3.567, -1.477]
_fcool(T) = 1.0 - 10^(sum([a_ * (log10(T))^(i - 1) for (i, a_) in enumerate(a)]))
function fcool(T) 
     if T > 10.0 
    return  1.0 - _fcool(10.0) 
elseif T < 0.0
    return 0.0
else 
    return 1.0 - _fcool(T)
end
end

function compute_qrad(p::StangebyModelParameters, r ; Te_init=10.0) 
    r.T_target = get_Te_target(p.sol.T_up, r.q_parallel_omp, p.target.L_para, p.integral.f_cond, p.integral.κ0; Te_init)
    return (1 - fcool(r.T_target)) * r.q_parallel_omp
end


function compute_zeff_up(par::StangebyModelParameters)
    zeff = ADAS.get_Zeff(par.sol.imp)
    return zeff(par.sol.f_imp, par.sol.n_up, par.sol.T_up)
end

"""
    summary(model::StangebyModel)

Print summary of StangebyModel setup and simulation results
"""
function Base.summary(model::StangebyModel)
    p = model.parameters
    r = model.results
    printfmtln("Upstream")
    printfmtln("├─ {:<22} = {:.2f} MW", "P_SOL", p.plasma.P_SOL / 1e6)
    printfmtln("├─ {:<22} = {:.2f} eV", "Te_up", p.sol.T_up)
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
    printfmtln("Stangeby model output")
    printfmtln("├─ {:<22} = {:.2f} MW/m^2", "q_poloidal_omp", r.q_poloidal_omp / 1e6)
    printfmtln("├─ {:<22} = {:.2f} MW/m^2", "q_parallel_omp", r.q_parallel_omp / 1e6)
    printfmtln("├─ {:<22} = {:.2f} MW/m^2 ", "q_rad", (r.q_rad) / 1e6)
    printfmtln("├─ {:<22} = {:.2f} MW/m^2 ", "q_rad_effective", (r.q_parallel_omp - r.q_parallel_target_unprojected) / 1e6)
    printfmtln("├─ {:<22} = {:.2f} MW/m^2", "q_parallel_target", r.q_parallel_target_unprojected / 1e6)
    printfmtln("├─ {:<22} = {:.2f} MW/m^2 ", "q_perp_target", r.q_perp_target / 1e6)
    printfmtln("└─ {:<22} = {:.2f} MW/m^2 ", "q_perp_target_spread", r.q_perp_target_spread / 1e6)
end
include("Stangeby_from_dd.jl")
end