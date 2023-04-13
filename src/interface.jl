



import IMASDD
mutable struct DivertorHeatFlux{D<:Union{Nothing,IMASDD.dd}, M<: DivertorHeatFluxModel}
    dd :: D
    model :: M
end
DivertorHeatFlux(; model::Symbol=:lengyel) = DivertorHeatFlux(nothing, model)
DivertorHeatFlux(dd::Union{Nothing,IMASDD.dd}; model::Symbol=:lengyel) = DivertorHeatFlux(dd, model)
DivertorHeatFlux(dd::Union{Nothing,IMASDD.dd}, model::Symbol) = DivertorHeatFlux(dd, DivertorHeatFluxModel(model))

function DivertorHeatFluxModel(model::Symbol)
    if model == :lengyel 
        return  LengyelModel()
    else
        error()
    end 
end
function DivertorHeatFluxModel(model::Symbol,setup)
    if model == :lengyel 
        return  LengyelModel(setup.lengyel)
    else
        error()
    end 
end

function DivertorHeatFlux(dd::IMASDD.dd,par::DivertorHeatFluxParameters) 
    return DivertorHeatFlux(dd,par.model,par.setup)
end
setup_model(dhf::DivertorHeatFlux) = setup_model(dhf.model,dhf.dd)



