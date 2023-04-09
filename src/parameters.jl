



using SimulationParameters: AbstractParameters, Entry, Switch
import .LengyelHeatFluxModel:LengyelModelParameters

Base.@kwdef mutable struct DivertorHeatFluxModelParameters{T} <: AbstractParameters where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    lengyel ::  LengyelModelParameters{T} = LengyelModelParameters()
end


abstract type FUSEparameters__ActorDivertorHeatFlux{T} <: AbstractParameters end
Base.@kwdef mutable struct DivertorHeatFluxParameters{T<:Real} <: FUSEparameters__ActorDivertorHeatFlux{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set 
    model:: Switch{Symbol} = Switch(Symbol, [:lengyel], "-", ""; default=:lengyel)
    setup :: DivertorHeatFluxModelParameters{T} = DivertorHeatFluxModelParameters()
end

FUSEparameters__ActorDivertorHeatFlux{T}(args...; kw...) where {T<:Real} = DivertorHeatFluxParameters{T}(args...; kw...)