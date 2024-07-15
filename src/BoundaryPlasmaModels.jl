module BoundaryPlasmaModels
using IMAS

abstract type DivertorHeatFluxModel end

include(joinpath("LengyelModel", "Lengyel_parameters.jl"))
include(joinpath("LengyelModel", "LengyelHeatFluxModel.jl"))
include(joinpath("LengyelModel", "Lengyel_from_dd.jl"))

include("parameters.jl")

include("interface.jl")

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__, all=false, imported=false) if name != Symbol(@__MODULE__)]

end
