"""
    setup_model(
        boundary_plasma_model::LengyelModel,
        eqt::IMAS.equilibrium__time_slice,
        cp1d::IMAS.core_profiles__profiles_1d,
        sol::Vector{IMAS.OpenFieldLine},
        P_SOL::Real,
        λ_omp::Real;
        strike_index::Int,
        imp::Vector{Symbol},
        f_imp::Vector{<:Real},
        f_spread_pfr::Real=1.0)

`strike_index` selects the strike point to USED

  - 0: automatic strike point selection favoring the outermost strike point
  - 1: strike point at the 1st location in the `sol` OpenFieldLine
  - anything else: strike point at the `end` location in the `sol` OpenFieldLine
"""
function setup_model(
    boundary_plasma_model::StangebyModel,
    target::IMAS.divertors__divertor___target,
    eqt::IMAS.equilibrium__time_slice,
    cp1d::IMAS.core_profiles__profiles_1d,
    sol1::IMAS.OpenFieldLine;
    impurities::Vector{Symbol},
    impurities_fraction::Vector{<:Real},
    heat_spread_factor::Real=1.0)

    @assert (length(impurities) == length(impurities_fraction))

    boundary_plasma_model.parameters.plasma.P_SOL = @ddtime(target.power_conducted.data) + @ddtime(target.power_convected.data)
    boundary_plasma_model.parameters.plasma.R_omp = sol1.r[sol1.midplane_index]
    boundary_plasma_model.parameters.plasma.Ip = eqt.global_quantities.ip
    boundary_plasma_model.parameters.plasma.κ = eqt.boundary.elongation
    boundary_plasma_model.parameters.plasma.ϵ = eqt.boundary.minor_radius / eqt.boundary.geometric_axis.r
    boundary_plasma_model.parameters.plasma.Bt_omp = sol1.Bt[sol1.midplane_index]
    boundary_plasma_model.parameters.plasma.Bpol_omp = sol1.Bp[sol1.midplane_index]

    boundary_plasma_model.parameters.sol.T_up = cp1d.electrons.temperature[end]
    boundary_plasma_model.parameters.sol.λ_omp = target.two_point_model[].sol_heat_decay_length

    boundary_plasma_model.parameters.target.f_omp2target_expansion = @ddtime(target.flux_expansion.data)
    λ_target = boundary_plasma_model.parameters.sol.λ_omp * boundary_plasma_model.parameters.target.f_omp2target_expansion
    boundary_plasma_model.parameters.target.R = @ddtime(target.wetted_area.data) / (λ_target * 2π) #why not strike point location????
    boundary_plasma_model.parameters.target.L_para = sol1.s[end]
    boundary_plasma_model.parameters.target.f_spread_pfr = heat_spread_factor
    boundary_plasma_model.parameters.target.α_sp = @ddtime(target.tilt_angle_tor.data)
    return boundary_plasma_model.parameters.target.θ_sp = @ddtime(target.tilt_angle_pol.data)
end