module BoundaryPlasmaModels
using IMAS
using SimulationParameters: AbstractParameters, Entry, Switch

abstract type SOLBoundaryModel end
abstract type SOLBoundaryParameters{T} <: AbstractParameters{T} end

# --- backward-compat: MUST be before includes ---
const DivertorHeatFluxModel = SOLBoundaryModel
const DivertorHeatFluxModelParameters{T} = SOLBoundaryParameters{T} where {T}

include("StangebyModel/StangebyHeatFluxModel.jl")
include("LengyelModel/LengyelHeatFluxModel.jl")
include("SLCoupledBoundary/SLCoupledModel.jl")   


using .StangebyHeatFluxModel
using .LengyelHeatFluxModel
using .SLCoupledFullModel 

# --- parameters --- #

Base.@kwdef mutable struct SOLBoundaryParametersShell{T<:Real} <: AbstractParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch{Symbol}([:slcoupled,:lengyel,:stangeby], "-", "SOL boundary model"; default=:slcoupled)
    setup::SOLBoundaryParameters{T} = SLCoupledParameters{T}()
end
const DivertorHeatFluxParameters = SOLBoundaryParametersShell   

const _setup_dic = Dict(:slcoupled => SLCoupledParameters, :lengyel => LengyelModelParameters, :stangeby => StangebyModelParameters)

function Base.setproperty!(par::SOLBoundaryParametersShell{T}, prop::Symbol, v) where {T}
    if prop == :model
        setproperty!(getfield(par,:model), :value, v)
        setfield!(par, :setup, _setup_dic[v]{T}())
        println("setup model parameters for model `$v`")
    elseoka
        setfield!(par, prop, v)
    end
end
# --- entry point --- #
function SOLBoundaryModel(model::Symbol)
    if model == :lengyel
        return LengyelModel()
    elseif model == :stangeby
        return StangebyModel()
    elseif model == :slcoupled
        return SLCoupledModel()
    else
        error("SOL boundary model `$model` is not recognized")
    end
end

SOLBoundaryModel(par::SOLBoundaryParametersShell) = SOLBoundaryModel(getfield(par, :model).value, par.setup)

function SOLBoundaryModel(model::Symbol, setup::SOLBoundaryParameters)
    if     model == :lengyel
        return LengyelModel(setup::LengyelModelParameters)
    elseif model == :stangeby
        return StangebyModel(setup::StangebyModelParameters)
    elseif model == :slcoupled
        return SLCoupledModel(setup::SLCoupledParameters)
    else
        error("SOL boundary model `$model` is not recognized")
    end
end

setup_model(model::LengyelHeatFluxModel.LengyelModel, args...; kw...) = LengyelHeatFluxModel.setup_model(model, args...; kw...)
setup_model(model::StangebyHeatFluxModel.StangebyModel, args...; kw...) = StangebyHeatFluxModel.setup_model(model, args...; kw...)
setup_model(model::SLCoupledFullModel.SLCoupledModel, args...; kw...) = SLCoupledFullModel.setup_model(model, args...; kw...)

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__, all=false, imported=false) if name != Symbol(@__MODULE__)]

end



