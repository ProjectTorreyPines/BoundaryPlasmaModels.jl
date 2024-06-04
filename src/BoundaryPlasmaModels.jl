module BoundaryPlasmaModels
using IMAS
using SimulationParameters: AbstractParameters, Entry, Switch

abstract type DivertorHeatFluxModel end
abstract type DivertorHeatFluxModelParameters{T} <: AbstractParameters{T} end
include("StangebyModel/StangebyHeatFluxModel.jl")
include("LengyelModel/LengyelHeatFluxModel.jl")

using .StangebyHeatFluxModel
using .LengyelHeatFluxModel

# --- parameters --- #

Base.@kwdef mutable struct DivertorHeatFluxParameters{T<:Real} <: AbstractParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch{Symbol}([:lengyel], "-", ""; default=:lengyel)
    setup::DivertorHeatFluxModelParameters{T} = LengyelModelParameters{T}()
end

# --- entry point --- #
function DivertorHeatFluxModel(model::Symbol)
    if model == :lengyel
        return LengyelModel()
    elseif model == :stangeby
        return StangebyModel()
    else
        error("Divertor heat flux model `$model` is not recognized")
    end
end

DivertorHeatFluxModel(par::DivertorHeatFluxParameters) = DivertorHeatFluxModel(par.model, par.setup)

function DivertorHeatFluxModel(model::Symbol, setup::DivertorHeatFluxModelParameters)
    if model == :lengyel
        return LengyelModel(setup)
    elseif model == :stangeby
        return StangebyModel(setup)
    else
        error()
    end
end



end
