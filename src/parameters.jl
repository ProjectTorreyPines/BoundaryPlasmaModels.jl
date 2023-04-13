



using SimulationParameters: AbstractParameters, Entry, Switch

Base.@kwdef mutable struct DivertorHeatFluxModelParameters{T} <: AbstractParameters where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    lengyel ::  LengyelModelParameters{T} = LengyelModelParameters{T}()
end


abstract type FUSEparameters__ActorDivertorHeatFluxTarget{T} <: AbstractParameters end
Base.@kwdef mutable struct DivertorHeatFluxParameters{T<:Real} <: FUSEparameters__ActorDivertorHeatFluxTarget{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set 
    model:: Switch{Symbol} = Switch(Symbol, [:lengyel], "-", ""; default=:lengyel)
    setup :: DivertorHeatFluxModelParameters{T} = DivertorHeatFluxModelParameters{T}()
end

FUSEparameters__ActorDivertorHeatFluxTarget{T}(args...; kw...) where {T<:Real} = DivertorHeatFluxParameters{T}(args...; kw...)