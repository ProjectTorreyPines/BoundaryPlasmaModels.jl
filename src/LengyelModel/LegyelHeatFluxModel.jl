using LaTeXStrings
import NumericalIntegration
import ADAS
import IMASDD
import FusionGeometryTools
using Formatting

include("parameters.jl")

mutable struct LengyelModelResults{T}
    q_poloidal_omp::T
    q_parallel_omp::T
    q_rad::T
    q_parallel_target_unprojected::T
    q_parallel_target_projected::T
    q_perp_target::T
    q_perp_target_spread::T
    zeff_up::T
end

function Base.show(io::IO, r::LengyelModelResults)
    for f in propertynames(r)
        println(io, "$f: $(getfield(r,f))")
    end
end

LengyelModelResults() = LengyelModelResults([0.0 for d in fieldnames(LengyelModelResults)]...)

mutable struct LengyelModel <: DivertorHeatFluxModel
    parameters::LengyelModelParameters
    results::LengyelModelResults
end

LengyelModel() = LengyelModel(LengyelModelParameters{Real}())

LengyelModel(par::LengyelModelParameters) = LengyelModel(par, LengyelModelResults())

(model::LengyelModel)() = model.results = compute_lengyel_model(model.parameters)

function setup_model(model::LengyelModel, dd::IMAS.dd)
    set_heatfluxmodel_plasma_parameters_from_dd(model, dd)
    set_heatfluxmodel_sol_parameters_from_dd(model, dd)
    set_λ_omp(model)
    set_heatfluxmodel_target_parameters_from_dd(model, dd)
end

function set_heatfluxmodel_sol_parameters(model::LengyelModel, dd::IMAS.dd)
    if model.parameters.sol.method == :dd
        set_heatfluxmodel_sol_parameters_from_dd(model, dd)
    elseif model.parameters.sol.method == :manual
    else
        error()
    end
end

function set_heatfluxmodel_plasma_parameters(model::LengyelModel, dd::Union{IMAS.dd,Nothing})
    if model.parameters.plasma.method == :dd
        set_heatfluxmodel_plasma_parameters_from_dd(model, dd)
    elseif model.parameters.plasma.method == :manual
    else
        error("model.parameters.plasma.method = $(model.parameters.plasma.method)")
    end
end

function set_heatfluxmodel_plasma_parameters_from_dd(model::LengyelModel, dd::IMAS.dd)
    d = Dict(
        :P_SOL => get_PSOL(dd),
        :R_omp => FusionGeometryTools.get_R_omp(dd),
        :Ip => get_Ip(dd),
        :κ => get_κ(dd),
        :ϵ => get_ϵ(dd),
        :Bt_omp => FusionGeometryTools.get_Bt_omp(dd),
        :Bpol_omp => FusionGeometryTools.get_Bpol_omp(dd)
    )
    for (k, v) in d
        setproperty!(model.parameters.plasma, k, v)
    end
end

function set_heatfluxmodel_sol_parameters_from_dd(model::LengyelModel, dd::IMAS.dd)
    d = Dict(
        :n_up => get_omp_density(dd),
        :T_up => get_omp_Te(dd)
    )
    p = model.parameters.sol
    for (k, v) in d
        setproperty!(p, k, v)
    end
end

function set_heatfluxmodel_target_parameters_from_dd(model::LengyelModel, dd::IMAS.dd)
    loc = Symbol.(split(string(model.parameters.target.location), "_"))
    d = Dict(
        :λ_target => get_projected_λ_omp_target(model.parameters.sol.λ_omp, dd; location=loc),
        :f_omp2target_expension => get_projected_λ_omp_target(model.parameters.sol.λ_omp, dd; location=loc) / model.parameters.sol.λ_omp,
        :α_sp => get_α_sp(dd, loc) * 180 / pi,
        :f_pol_projection => tan(get_α_sp(dd, loc)),
        :θ_sp => get_θ_sp(dd, loc) * 180 / pi,
        :f_perp_projection => 1 / cos(get_θ_sp(dd, loc))
    )
    p = model.parameters.target
    for (k, v) in d
        setproperty!(p, k, v)
    end
end

function set_λ_omp(model)
    if model.parameters.sol.λ_omp_scaling isa Number
        model.parameters.sol.λ_omp = model.parameters.sol.λ_omp_scaling
    else
        model.parameters.sol.λ_omp = get_λ_omp(model.parameters)
    end
end

function get_λ_omp(parameters)
    if parameters.sol.λ_omp_scaling == :eich
        return λq_eich(parameters.plasma)
    elseif parameters.sol.λ_omp_scaling == :none
        return parameters.sol.λ_omp
    else
        error()
    end
end

""" Perform weighted cooling rate integral over specified temperature interval
# Inputs:  Tmin   minimum temperature for integral (eV)
#          Tmax   maximum temperature for integral (eV)
#          Texp   Exponent for temperature weighting
#          Lexp   Exponent for cooling rate weighting
#          Zimp   Z of impurity (Ar: 18; Kr: 36; Xe: 54)
"""
function V_legyel_ADAS(Tmin::Float64, Tmax::Float64, f_imp::Float64, imp::Union{String,Symbol}; N::Int64=500, ne::Float64=1e20, Zeff_exp::Float64=-0.3, Texp::Float64=0.5, Lexp::Float64=1.0, κ0=2390.0)
    data = ADAS.get_cooling_rates(imp)
    zeff = ADAS.get_effective_charge(imp)
    Lz = data.Lztot
    T = collect(LinRange(Tmin, Tmax, N))
    int = [T_ .^ Texp .* zeff(f_imp, ne, T_) .^ (Zeff_exp) .* Lz(ne, T_) .^ Lexp for T_ in T]
    return sqrt.(NumericalIntegration.integrate(T, int) * f_imp * κ0 * 2)
end

function V_legyel_ADAS(Tmin::Float64, Tmax::Float64, f_imps::Vector{Float64}, imps::Vector{<:Union{String,Symbol}}; ne::Float64=1e20, Zeff_exp::Float64=-0.3, Texp::Float64=0.5, Lexp::Float64=1.0, κ0=2390.0, N::Int64=500)
    zeff = ADAS.get_effective_charge(imps)
    T = collect(LinRange(Tmin, Tmax, N))
    int = 0.0
    for (f_imp, imp) in zip(f_imps, imps)
        data = ADAS.get_cooling_rates(imp)
        int += sqrt.(NumericalIntegration.integrate(T, [T_ .^ Texp .* zeff(f_imps, ne, T_) .^ (Zeff_exp) .* data.Lztot(ne, T_) .^ Lexp for T_ in T]) * f_imp * κ0 * 2.0)
    end
    return int
end

V_legyel_ADAS(s, i) = V_legyel_ADAS(i.T_down, s.T_up, s.f_imp, s.imp; ne=s.n_up, Zeff_exp=i.Zeff_exp, Texp=i.Texp, Lexp=i.Lexp, κ0=i.κ0)

function compute_lengyel_model(par::LengyelModelParameters)
    r = LengyelModelResults()
    r.q_poloidal_omp = compute_q_poloidal_omp(par)
    r.q_parallel_omp = compute_q_parallel_omp(par)
    r.q_rad = compute_qrad(par)
    r.q_parallel_target_unprojected = sqrt(r.q_parallel_omp^2.0 - r.q_rad^2.0)
    r.q_parallel_target_projected = r.q_parallel_target_unprojected / par.target.f_omp2target_expension
    r.q_perp_target = r.q_parallel_target_projected * par.target.f_pol_projection
    r.q_perp_target_spread = r.q_perp_target * par.target.f_spread_pfr
    r.zeff_up = compute_zeff_up(par)
    return r
end

compute_q_parallel_omp(p::LengyelModelParameters) = compute_q_parallel_omp(p.plasma.P_SOL, p.plasma.R_omp, p.sol.λ_omp, p.plasma.Bpol_omp, p.plasma.Bt_omp)
compute_q_poloidal_omp(p::LengyelModelParameters) = compute_q_poloidal_omp(p.plasma.P_SOL, p.plasma.R_omp, p.sol.λ_omp)
compute_q_poloidal_omp(P_SOL::T, R::T, λ_q::T) where {T<:Float64} = begin
    @assert R > 0.0 && λ_q > 0
    P_SOL / (2.0 * pi * R * λ_q)
end
compute_q_parallel_omp(P_SOL::T, R::T, λ_q::T, Bpol::T, Bt::T) where {T<:Float64} = begin
    @assert (R > 0.0 && Bpol > 0 && Bt > 0.0 && λ_q > 0)
    P_SOL / (2.0 * pi * R * λ_q) / sin(atan(Bpol / Bt))
end

compute_qrad(p::LengyelModelParameters) = compute_qrad(p.sol, p.integral)
compute_qrad(s, i) = s.f_adhoc * s.n_up * s.T_up * V_legyel_ADAS(s, i)

function compute_zeff_up(par::LengyelModelParameters)
    zeff = ADAS.get_effective_charge(par.sol.imp)
    return zeff(par.sol.f_imp, par.sol.n_up, par.sol.T_up)
end

function show_summary(model::LengyelModel)
    printfmtln("---- Lengyel model for divertor: {}", model.parameters.target.location)
    printfmtln("└─ upstream parameters ")
    printfmtln("  └─ {:<10} = {:.2f} MW", "P_SOL", model.parameters.plasma.P_SOL / 1e6)
    printfmtln("  └─ {:<10} = {:.2f} eV", "Te_up", model.parameters.sol.T_up)
    printfmtln("  └─ {:<10} = {:.2e} m^-3", "ne_up", model.parameters.sol.n_up)
    printfmtln("  └─ {:<10} = {:.1f} T", "Bp_omp", model.parameters.plasma.Bpol_omp)
    printfmtln("  └─ {:<10} = {:.1f} m", "R_omp", model.parameters.plasma.R_omp)
    printfmtln("  └─ {:<10} = {:.1f} MW.T/m", "PBp/R", model.parameters.plasma.P_SOL / 1e6 / model.parameters.plasma.R_omp * model.parameters.plasma.Bpol_omp)
    printfmtln("  └─ {:<10} = {:.4f} m [scaling law:{}]", "λ_omp", model.parameters.sol.λ_omp, string(model.parameters.sol.λ_omp_scaling))
    printfmtln("└─ target parameters ")
    printfmtln("  └─ {:<10} = {:.4f} m ", "λ_target", model.parameters.target.λ_target)
    printfmtln("  └─ {:<10} = {:.1f} deg", "α_pitch", model.parameters.target.α_sp)
    printfmtln("  └─ {:<10} = {:.1f} deg", "θ_target", model.parameters.target.θ_sp)
    printfmtln("  └─ {:<10} = {:.1f} deg", "θ_target", model.parameters.target.θ_sp)
    printfmtln("  └─ {:<10} = {:.1f} m", "spread_pfr", model.parameters.target.f_spread_pfr)
    printfmtln("└─ impurities ")
    printfmtln("  └─ name | fraction ")
    for (i, f) in zip(model.parameters.sol.imp, model.parameters.sol.f_imp)
        printfmtln("  └─ {:<4} | {:<10.3f}  ", string(i), f)

    end
    printfmtln("  └─ {:<15} = {:.2f} ", "Zeff_up", model.results.zeff_up)
    printfmtln("└─ model output ")
    printfmtln("  └─ {:<15} = {:.2f} MW/m^2", "q_parallel_omp", model.results.q_parallel_omp / 1e6)
    printfmtln("  └─ {:<15} = {:.2f} MW/m^2", "q_poloidal_omp", model.results.q_poloidal_omp / 1e6)
    printfmtln("  └─ {:<15} = {:.2f} MW/m^2 ", "q_rad", (model.results.q_rad) / 1e6)
    printfmtln("  └─ {:<15} = {:.2f} MW/m^2 ", "q_rad_effective", (model.results.q_parallel_omp - model.results.q_parallel_target_unprojected) / 1e6)
    printfmtln("  └─ {:<15} = {:.2f} MW/m^2", "q_parallel_target", model.results.q_parallel_target_unprojected / 1e6)
    printfmtln("  └─ {:<15} = {:.2f} MW/m^2 ", "q_perp_target", model.results.q_perp_target / 1e6)
    printfmtln("  └─ {:<15} = {:.2f} MW/m^2 ", "q_perp_target_spread", model.results.q_perp_target_spread / 1e6)
end
