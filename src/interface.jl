



import IMASDD
mutable struct DivertorHeatFlux
    dd :: Union{Nothing,IMASDD.dd}
    model :: DivertorHeatFluxModel
end

DivertorHeatFluxModel(dd::IMASDD.dd; model::Symbol=:lengyel) = DivertorHeatFlux(dd, model)
DivertorHeatFluxModel(; model::Symbol=:lengyel) = DivertorHeatFlux(nothing, model)

DivertorHeatFlux(dd::IMASDD.dd, model::Symbol=:lengyel) = DivertorHeatFlux(dd, DivertorHeatFluxModel(model))


function DivertorHeatFluxModel(model::Symbol)
    if model == :lengyel 
        return  LengyelHeatFluxModel.LengyelModel()
    else
        error()
    end 
end

function DivertorHeatFlux(dd::IMASDD.dd,par::DivertorHeatFluxParameters) 
    # dhf DivertorHeatFlux(dd,par.model)
    dhf.model(par.setup)
    return dhf
end

setup_model(dhf::DivertorHeatFlux, setup) = setup_model(dhf.model,dhf.dd,setup)



