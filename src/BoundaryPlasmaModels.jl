module BoundaryPlasmaModels
using DocStringExtensions

abstract type DivertorHeatFluxModel end 

include("LengyelModel/LegyelHeatFluxModel.jl")
include("parameters.jl")
include("interface.jl")
include("basic.jl")

end
