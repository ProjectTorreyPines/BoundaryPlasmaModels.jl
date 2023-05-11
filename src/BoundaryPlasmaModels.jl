module BoundaryPlasmaModels
using DocStringExtensions
using RecipesBase
using IMAS

abstract type DivertorHeatFluxModel end 

include(joinpath("LengyelModel","Lengyel_parameters.jl"))
include(joinpath("LengyelModel","LengyelHeatFluxModel.jl"))
include(joinpath("LengyelModel","Lengyel_from_dd.jl"))

include("parameters.jl")
include("interface.jl")
include("scaling_laws.jl")
include("utils.jl")

end
