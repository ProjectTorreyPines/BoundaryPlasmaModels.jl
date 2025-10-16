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
    model::Switch{Symbol} = Switch{Symbol}([:lengyel,:stangeby], "-", ""; default=:lengyel)
    setup::DivertorHeatFluxModelParameters{T} = LengyelModelParameters{T}()
end

setup_model_dic = Dict(:lengyel => LengyelModelParameters, :stangeby => StangebyModelParameters)

function Base.setproperty!(par::DivertorHeatFluxParameters{T},prop::Symbol, v) where T
    if prop == :model
        setproperty!(getfield(par,:model),:value,v) 
        setfield!(par,:setup,setup_model_dic[v]{T}())
        println("setup model parameters for model `$v`") 
    else
        setfield!(par,prop,v)
    end
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
        error("Divertor heat flux model `$model` is not recognized")
    end
end

setup_model(model::LengyelModel, args...; kw...) = LengyelHeatFluxModel.setup_model(model, args...; kw...)
setup_model(model::StangebyModel, args...; kw...) = StangebyHeatFluxModel.setup_model(model, args...; kw...)


const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__, all=false, imported=false) if name != Symbol(@__MODULE__)]

end
