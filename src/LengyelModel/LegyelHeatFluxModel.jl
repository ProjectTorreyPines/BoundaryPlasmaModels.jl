
module LengyelHeatFluxModel
using LaTeXStrings
using RecipesBase
import NumericalIntegration
import ADAS
import IMASDD
import ..DivertorHeatFluxModel
import FusionGeometryTools
include("parameters.jl")

mutable struct LengyelModelResults{T}
    q_poloildal_omp :: T
    q_parallel_omp :: T 
    q_rad :: T 
    q_parallel_target :: T 
    q_parallel_target_expended::T 
    q_parallel_target_expended_spread :: T
    q_pol_target :: T
    q_perp_target :: T
end

LengyelModelResults() = LengyelModelResults([0.0 for d in fieldnames(LengyelModelResults)]...)

mutable struct LengyelModel <: DivertorHeatFluxModel
    parameters :: LengyelModelParameters
    results :: LengyelModelResults
end

LengyelModel() = LengyelModel(LengyelModelParameters{Real}(),LengyelModelResults())

(model::LengyelModel)() = model.results = compute_lengyel_model(model.parameters)


function setup_plasma_model(model::LengyelModel, dd::Union{IMASDD.dd, Nothing}, setup)
    set_heatfluxmodel_plasma_parameters(model,dd, setup.plasma)
    set_heatfluxmodel_sol_parameters(model,dd, setup.sol)
end

function set_heatfluxmodel_sol_parameters(model::LengyelModel, dd::Union{IMASDD.dd, Nothing}, setup_sol)
    if setup_sol.method == :auto
        set_heatfluxmodel_sol_parameters_from_dd(model,dd)
    elseif setup_sol.method ==:manual
    else
        error()
    end
    set_λ_omp_scaling(model,setup_sol)
    set_λ_omp(model)
end


get_P_SOL(a) = 0.0

function set_heatfluxmodel_plasma_parameters_from_dd(model::LengyelModel, dd:: IMASDD.dd)
d = Dict( 
    :P_SOL => get_P_SOL(dd),
    :R_omp => FusionGeometryTools.get_R_omp(dd.eqt),
    :Ip => FusionGeometryTools.get_ip(dd.eqt),
    :κ => FusionGeometryTools.get_κ(dd.eqt),
    :ϵ => FusionGeometryTools.get_ϵ(dd.eqt),
    :Bt_omp => FusionGeometryTools.get_Bt_omp(dd.eqt),
    :Bpol_omp => FusionGeometryTools.get_Bpol_omp(dd.eqt)
    )
    for f in propertynames(model.parameters.plasma)
        p = model.parameters.plasma
        setfield!(p,f,d[f])
    end
end

set_λ_omp_scaling(model,setup_sol) = model.parameters.sol.λ_omp_scaling = setup_sol.λ_omp_scaling 
function set_λ_omp(model) 
    if model.parameters.sol.λ_omp_scaling isa Number
        model.parameters.plasma.λ_omp = model.parameters.sol.λ_omp_scaling 
    else
        model.parameters.plasma.λ_omp = model.parameters.sol.λ_omp_scaling(model.parameters)
    end
end

function set_heatfluxmodel_sol_parameters_from_dd(model::LengyelModel, dd:: IMASDD.dd)
    d=Dict( 
        :Ndiv => get_ndiv(dd),
        :n_up => get_omp_density(dd),
        :T_up => get_omp_Te(dd)
        )
        p = model.parameters.sol
    for (k,v) in d 
        setfield!(p,k, v)
    end
end





""" Perform weighted cooling rate integral over specified temperature interval
# Inputs:  Tmin   minimum temperature for integral (eV)
#          Tmax   maximum temperature for integral (eV)
#          Texp   Exponent for temperature weighting
#          Lexp   Exponent for cooling rate weighting
#          Zimp   Z of impurity (Ar: 18; Kr: 36; Xe: 54)
"""
function V_legyel_ADAS(Tmin::Float64, Tmax::Float64, f_imp :: Float64, imp::Union{String,Symbol} ; N::Int64= 500, ne::Float64 = 1e20, Zeff_exp::Float64=-0.3,Texp::Float64=0.5,Lexp::Float64 =1.0,κ0 = 2390.0)
    data = ADAS.get_cooling_rates(imp);
    zeff = ADAS.get_effective_charge(imp);
    Lz = data.Lztot 
    T = collect(LinRange(Tmin,Tmax,N))
    int = [T_ .^ Texp .* zeff(f_imp,ne,T_).^(Zeff_exp) .* Lz(ne,T_).^ Lexp for T_ in T]
    return sqrt.(NumericalIntegration.integrate(T,int) * f_imp * κ0)
end

function V_legyel_ADAS(Tmin::Float64, Tmax::Float64, f_imp :: Vector{Float64}, imp::Vector{<:Union{String,Symbol}} ; ne::Float64 = 1e20, Zeff_exp::Float64=-0.3,Texp::Float64=0.5,Lexp::Float64 =1.0,κ0 = 2390.0,N::Int64= 500)
    
    zeff = ADAS.get_effective_charge(imp);
    T = collect(LinRange(Tmin,Tmax,N))
    int = 0.0
    for (f,i) in zip(f_imp,i)
        data = ADAS.get_cooling_rates(i); 
        int += sqrt.(NumericalIntegration.integrate(T,[T_ .^ Texp .* zeff(f_imp,ne,T_).^(Zeff_exp) .* data.Lztot(ne,T_).^ Lexp for T_ in T])* f * κ0)
    end
    return int
end
V_legyel_ADAS(s,i) = V_legyel_ADAS(i.T_down,s.T_up,s.f_imp, s.imp; ne = s.n_up, Zeff_exp=i.Zeff_exp,Texp=i.Texp,Lexp=i.Lexp,κ0=i.κ0)
function compute_lengyel_model(par::LengyelModelParameters)
    r = LengyelModelResults()
    r.q_poloildal_omp = compute_q_poloidal_omp(par)
    r.q_parallel_omp = compute_q_parallel_omp(par)
    r.q_rad = compute_qrad(par) 
    r.q_parallel_target = sqrt(r.q_parallel_omp^2.0 - r.q_rad^2.0)
    r.q_parallel_target_expended = r.q_parallel_target_omp * par.target.f_omp2target_expension 
    r.q_parallel_target_expended_spread = r.q_parallel_target_expended * f_spread_pfr
    r.q_pol_target = r.q_parallel_target_projected * par.target.f_pol_projection
    r.q_perp_target = r.q_pol_target  * par.target.f_perp_projection
    return r
end




compute_q_parallel_omp(p::LengyelModelParameters) = compute_q_parallel_omp(p.plasma.P_SOL, p.plasma.R_omp,p.sol.λ_omp, p.plasma.Bpol_omp, p.plasma.Bt_omp)
compute_q_poloidal_omp(p::LengyelModelParameters) = compute_q_poloidal_omp(p.plasma.P_SOL, p.plasma.R_omp,p.sol.λ_omp)
compute_q_poloidal_omp(P_SOL::T, R::T,λ_q::T) where {T<: Float64}= P_SOL/(2.0*pi*R* λ_q)
compute_q_parallel_omp(P_SOL::T, R::T,λ_q::T, Bpol::T, Bt::T) where {T<: Float64}= P_SOL/(2.0*pi*R* λ_q)/sin(atan(Bpol/Bt))
compute_qrad(p::LengyelModelParameters) = compute_qrad(p.sol,p.integrale)
compute_qrad(s,i) = s.f_adhoc * s.n_up * s.T_up * V_legyel_ADAS(s,i)

# μ₀ = 4.e-7 * pi
# compute_Tup( q_parallel::T, ℒ_parallel::T,κ₀::T) = (7.*q_parallel*ℒ_parallel/(2.*κ₀))^(2./7.) 
# compute_q_star(plasma::LengyelModel)
# compute_q_star(R::T,κ::T,Ip::T,ϵ::T,Bt::T) where {T<:Real} = pi*(ϵ*R)^2. * (1. + κ^2.) * Bt/(μ₀*R*Ip*1.e6)
# compute_ℒ_omp_target(plasma :: GASC0DHeatFluxPlasmaParameters) = compute_ℒ_omp_target(plasma.R,plasma.κ,plasma.Ip,plasma.ϵ,plasma.Bt)
# compute_ℒ_omp_target(R::T,κ::T,Ip::T,ϵ::T,Bt::T) where {T<:Real} = 4.33 * compute_q⋆(R,κ,Ip,ϵ,Bt) * R
compute_Aperp(R::T,λ_q::T, Bpol::T, Bt::T) where {T<:Real} = 2.0*pi*R* λ_q 

# compute_κ₀(Zeff) where {T<:Real} = 2390. * Zeff^(-0.3) #mks units with T in eV (Stacey Fusion Plasma Physics pg. 379)  - Using 14. for ln(gamma)
# compute_q_rad_ADAS(n_up::T, T_up::T, imp::GASC0DHeatFluxSOLImpurityParameters) where {T<:Real} =    n_up * T_up *  sqrt(2. * compute_κ₀(imp.Zeff) *(imp.fzHe*Lint_ADAS(0.,T_up,0.5,1.0,imp.He)+
#                                        imp.fz1Sol*Lint_ADAS(0.,T_up,0.5,1.0,imp.imp1)+
#                                        imp.fz2Sol*Lint_ADAS(0.,T_up,0.5,1.0,imp.imp2) ))




# Lengyel_model(par::LengyelModelParameters) = Lengyel_model(par.plasma,par.sol,par.imp)
# function Lengyel_model(plasma, sol::LengyelModelSOLParameters)
#     ℒ_parallel = compute_ℒ_omp2target(plasma)
#     q_parallel_omp = compute_q_parallel_omp(plasma)
#     κ₀ = compute_κ₀(imp.Zeff)
#     T_up = compute_Tup( q_parallel_omp, ℒ_parallel,κ₀)
#     q_rad = compute_q_rad_ADAS(n_up, T_up, imp) 
#     q_target = sqrt(q_parallel_omp^2. - q_rad^2.)   #Expected heat flux without geometric considerations
#     P_rad = q_rad * A_perp * Ndiv * 1.e-6
#     return LengyelModelResults(q_parallel, ℒ_parallel, q_rad, q_target, T_up, P_rad) #, lambda_int)
# end




# function compute_fz_Reinke(P_SOL, R, Bt, Bpol, Ip, ϵ, κ , n_up, Zimp):

#     λ_q = 1.35 * P_SOL^(-0.02) *R ^0.04* Bpol^(-0.92) * ϵ^0.42 * 1.e-3
#     q_parallel = P_SOL*1.e6 * sqrt(Bt^2.+Bpol^2.) / (2.*pi*R*λ_q*Bpol)
#     #q_parallel2 = Psol*1.e6/(2.*np.pi*R*lambdaq*np.sin(np.arctan(Bpol/Bt)))

#     q_start = compute_q_star(R,κ,Ip,ϵ,Bt)
#     ℒ_parallel = 4.33*q_star*R

#     #q_parallel2 = 117.9*Psol*np.sqrt(Bt**2.+Bpol**2.)/ (R * epsilon**0.42)*1.e6

#     κ₀ = compute_κ₀(Zeff)
#     #kappa_0 = 3.09e3 / Zimp /14.  
#     T_up = compute_T_up( q_∥∥, ℒ_∥∥,κ₀)    #eV

#     fz = q_parallel^2./(2.*κ₀*n_up^2.*T_up^2.*Lint(0.,T_up,0.5,1.0,Zimp))

#     return fz
# end

end