function DivertorHeatFluxModel(model::Symbol)
    if model == :lengyel
        return LengyelModel()
    else
        error("Divertor heat flux model `$model` is not recognized")
    end
end

function DivertorHeatFluxModel(par::DivertorHeatFluxParameters)
    return DivertorHeatFluxModel(par.model, par.setup)
end

function DivertorHeatFluxModel(model::Symbol, setup::DivertorHeatFluxModelParameters)
    if model == :lengyel
        return LengyelModel(setup.lengyel)
    else
        error()
    end
end
