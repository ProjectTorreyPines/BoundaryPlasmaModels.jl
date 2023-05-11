function setup_model(model::LengyelModel, dd::IMAS.dd)
    set_heatfluxmodel_plasma_parameters_from_dd(model, dd)
    set_heatfluxmodel_sol_parameters_from_dd(model, dd)
    set_λ_omp_from_dd(model, dd)
    set_heatfluxmodel_target_parameters_from_dd(model, dd)
end

function set_heatfluxmodel_plasma_parameters_from_dd(model::LengyelModel, dd::IMAS.dd)
    d = Dict(
        :P_SOL => get_PSOL(dd),
        :R_omp => FusionGeometryTools.get_R_omp(dd),
        :Ip => get_Ip(dd),
        :κ => get_κ(dd),
        :ϵ => get_ϵ(dd),
        :Bt_omp => FusionGeometryTools.get_Bt_omp(dd),
        :Bpol_omp => FusionGeometryTools.get_Bpol_omp(dd)
    )
    for (k, v) in d
        setproperty!(model.parameters.plasma, k, v)
    end
end

function set_heatfluxmodel_sol_parameters_from_dd(model::LengyelModel, dd::IMAS.dd)
    d = Dict(
        :n_up => get_omp_density(dd),
        :T_up => get_omp_Te(dd)
    )
    p = model.parameters.sol
    for (k, v) in d
        setproperty!(p, k, v)
    end
end

function set_heatfluxmodel_target_parameters_from_dd(model::LengyelModel, dd::IMAS.dd)
    loc = Symbol.(split(string(model.parameters.target.location), "_"))
    d = Dict(
        :λ_target => get_projected_λ_omp_target(model.parameters.sol.λ_omp, dd; location=loc),
        :f_omp2target_expension => get_projected_λ_omp_target(model.parameters.sol.λ_omp, dd; location=loc) / model.parameters.sol.λ_omp,
        :α_sp => get_α_sp(dd, loc) * 180 / pi,
        :f_pol_projection => tan(get_α_sp(dd, loc)),
        :θ_sp => get_θ_sp(dd, loc) * 180 / pi,
        :f_perp_projection => 1 / cos(get_θ_sp(dd, loc))
    )
    p = model.parameters.target
    for (k, v) in d
        setproperty!(p, k, v)
    end
end

function set_λ_omp_from_dd(model, dd)
    model.parameters.sol.λ_omp = IMAS.widthSOL_eich(dd)
end
