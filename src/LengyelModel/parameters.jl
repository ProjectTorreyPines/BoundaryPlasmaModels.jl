
using SimulationParameters: AbstractParameters, Entry, Switch

Base.@kwdef mutable struct LengyelModelPlasmaParameters{T} <: AbstractParameters where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    method::Switch{Symbol} = Switch(Symbol, [:dd,:manual,:dict], "-", "Method to define the baffle design point"; default=:manual)
    P_SOL :: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=0.)
    R_omp :: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=0.)
    Ip :: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=0.)
    κ :: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=0.)
    ϵ :: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=0.)
    Bt_omp :: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=0.)
    Bpol_omp:: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=0.)
end

Base.@kwdef mutable struct LengyelModelSOLParameters{T} <: AbstractParameters where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    method::Switch{Symbol} = Switch(Symbol, [:dd,:manual,:dict], "-", "Method to define the baffle design point"; default=:manual)
    λ_omp_scaling::Switch{Union{Symbol,Float64}} = Switch(Union{Symbol,Float64}, [:eich,:none], "-", "Method to define the baffle design point"; default=:eich)
    n_up :: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=0.)
    T_up ::Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=0.)
    λ_omp ::Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=0.003)
    f_imp :: Entry{Vector{Float64}} = Entry(Vector{Float64}, "m", "Buffer between plasma and wall"; default=Vector{Float64}())
    imp :: Entry{Vector{Symbol}} = Entry(Vector{Symbol}, "m", "Buffer between plasma and wall"; default=Vector{Symbol}())
    tau_imp :: Entry{Vector{Float64}} = Entry(Vector{Float64}, "m", "Buffer between plasma and wall"; default=Vector{Float64}())
    f_adhoc :: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=1.0)
end

Base.@kwdef mutable struct LengyelModelTargetParameters{T} <: AbstractParameters where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    f_omp2target_expension :: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=1.0)
    f_perp_projection :: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=1.0)
    f_pol_projection :: Entry{T} = Entry(T, "m", "poloidal projection factor"; default=1.0)
    α_sp :: Entry{T} = Entry(T, "deg", "Pitch angle at the outer strike point"; default=3.0)
    θ_sp :: Entry{T} = Entry(T, "deg", "Poloidal angle of the target at the outer strike point"; default=90.0)
    f_spread_pfr :: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=1.0)
    λ_target ::Entry{T} = Entry(T, "m", "heat flux width at the target"; default=0.000)
end


Base.@kwdef mutable struct LengyelIntegralParameters{T} <: AbstractParameters where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    method :: Switch{Symbol} = Switch(Symbol, [:auto,:manual], "-", "Method to define the baffle design point"; default=:manual)
    T_down :: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=0.0)
    Zeff_exp :: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=-0.3)
    Texp :: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=0.5)
    Lexp :: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=1.0)
    κ0 :: Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=2390.0)
end

# Base.@kwdef mutable struct DivertorHeatFluxPlasmaParameters{T} <: AbstractParameters where {T<:Real}
#     _parent::WeakRef = WeakRef(nothing)
#     _name::Symbol = :not_set

# end

# Base.@kwdef mutable struct DivertorHeatFluxSOLParameters{T} <: AbstractParameters where {T<:Real}
#     _parent::WeakRef = WeakRef(nothing)
#     _name::Symbol = :not_set
    
# end
Base.@kwdef mutable struct LengyelModelParameters{T} <: AbstractParameters where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    sol :: LengyelModelSOLParameters{T} = LengyelModelSOLParameters{T}()
    plasma :: LengyelModelPlasmaParameters{T} = LengyelModelPlasmaParameters{T}()
    target :: LengyelModelTargetParameters{T} = LengyelModelTargetParameters{T}()
    integrale :: LengyelIntegralParameters{T} = LengyelIntegralParameters{T}()
end
