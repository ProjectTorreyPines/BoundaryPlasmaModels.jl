module BoundaryPlasmaModels
using DocStringExtensions
using RecipesBase
using IMAS

abstract type DivertorHeatFluxModel end 

include("LengyelModel/LegyelHeatFluxModel.jl")
include("parameters.jl")
include("interface.jl")
include("scaling_laws.jl")
include("utils.jl")

end
