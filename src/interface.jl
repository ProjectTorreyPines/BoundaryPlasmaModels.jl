"""
    DivertorHeatFluxModel(model::Symbol)

Pick a model. Supported [:lengyel]
"""
function DivertorHeatFluxModel(model::Symbol)
    if model == :lengyel
        return LengyelModel()
    else
        error("Divertor heat flux model `$model` is not recognized")
    end
end

"""
    DivertorHeatFluxModel(par::DivertorHeatFluxParameters)

Run a DivertorHeatFluxModel starting from a DivertorHeatFluxParameters
"""
function DivertorHeatFluxModel(par::DivertorHeatFluxParameters)
    return DivertorHeatFluxModel(par.model, par.setup)
end

"""
    DivertorHeatFluxModel(model::Symbol, setup::DivertorHeatFluxModelParameters)

Run a DivertorHeatFluxModel starting from a model and setup
"""
function DivertorHeatFluxModel(model::Symbol, setup::DivertorHeatFluxModelParameters)
    if model == :lengyel
        return LengyelModel(setup.lengyel)
    else
        error("Divertor heat flux model `$model` is not recognized")
    end
end

export DivertorHeatFluxModel