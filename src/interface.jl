# mutable struct DivertorHeatFlux{D<:Union{Nothing,IMAS.dd}, M<: DivertorHeatFluxModel}
#     model :: M
# end
DivertorHeatFlux(; model::Symbol=:lengyel) = DivertorHeatFlux(nothing, model)

DivertorHeatFlux(dd::Union{Nothing,IMAS.dd}; model::Symbol=:lengyel) = DivertorHeatFlux(dd, model)

DivertorHeatFlux(dd::Union{Nothing,IMAS.dd}, model::Symbol) = DivertorHeatFlux(dd, DivertorHeatFluxModel(model))

function DivertorHeatFluxModel(model::Symbol)
    if model == :lengyel
        return LengyelModel()
    else
        error()
    end
end

function DivertorHeatFluxModel(model::Symbol, setup::DivertorHeatFluxModelParameters)
    if model == :lengyel
        return LengyelModel(setup.lengyel)
    else
        error()
    end
end

DivertorHeatFluxModel(par::DivertorHeatFluxParameters) = DivertorHeatFluxModel(par.model, par.setup)

function export2dd(dd::IMAS.dd, model::DivertorHeatFluxModel; kw...)
end