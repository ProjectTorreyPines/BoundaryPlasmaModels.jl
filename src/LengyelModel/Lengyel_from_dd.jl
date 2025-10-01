
"""
    setup_model(
        boundary_plasma_model::LengyelModel,
        λ_omp::Real,
        eqt::IMAS.equilibrium__time_slice,
        cp1d::IMAS.core_profiles__profiles_1d,
        cs::IMAS.core_sources,
        sol1::IMAS.OpenFieldLine,
        impurities::Vector{Symbol},
        impurities_fraction::Vector{<:Real},
        heat_spread_factor::Real=1.0)

`strike_index` selects the strike point to USED
* 0: automatic strike point selection favoring the outermost strike point
* 1: strike point at the 1st location in the `sol` OpenFieldLine
* anything else: strike point at the `end` location in the `sol` OpenFieldLine
"""
function setup_model(
    boundary_plasma_model::LengyelModel,
    λ_omp::Real,
    eqt::IMAS.equilibrium__time_slice,
    cp1d::IMAS.core_profiles__profiles_1d,
    cs::IMAS.core_sources,
    sol1::IMAS.OpenFieldLine;
    impurities::Vector{Symbol},
    impurities_fraction::Vector{<:Real},
    heat_spread_factor::Real=1.0)

    @assert (length(impurities) == length(impurities_fraction))

    n_targets = 2
    boundary_plasma_model.parameters.plasma.P_SOL = IMAS.power_sol(cs, cp1d) / n_targets
    boundary_plasma_model.parameters.plasma.R_omp = sol1.r[sol1.midplane_index]
    boundary_plasma_model.parameters.plasma.Ip = eqt.global_quantities.ip
    boundary_plasma_model.parameters.plasma.κ = eqt.boundary.elongation
    boundary_plasma_model.parameters.plasma.ϵ = eqt.boundary.minor_radius / eqt.boundary.geometric_axis.r
    boundary_plasma_model.parameters.plasma.Bt_omp = sol1.Bt[sol1.midplane_index]
    boundary_plasma_model.parameters.plasma.Bpol_omp = sol1.Bp[sol1.midplane_index]

    boundary_plasma_model.parameters.sol.n_up = cp1d.electrons.density_thermal[end]
    boundary_plasma_model.parameters.sol.T_up = cp1d.electrons.temperature[end]
    boundary_plasma_model.parameters.sol.λ_omp = λ_omp
    boundary_plasma_model.parameters.sol.imp = impurities
    boundary_plasma_model.parameters.sol.f_imp = impurities_fraction

    boundary_plasma_model.parameters.target.f_omp2target_expansion = sol1.total_flux_expansion[strike_index]
    λ_target = boundary_plasma_model.parameters.sol.λ_omp * boundary_plasma_model.parameters.target.f_omp2target_expansion
    boundary_plasma_model.parameters.target.R = sol1.r[strike_index]
    boundary_plasma_model.parameters.target.f_spread_pfr = heat_spread_factor
    boundary_plasma_model.parameters.target.α_sp = atan(sol1.Bp[strike_index] / sol1.Bt[strike_index])
    boundary_plasma_model.parameters.target.θ_sp = sol1.strike_angles[strike_index == 1 ? 1 : 2]
end


export setup_model

# questions:
# NO R_target ?
# Heat flux treated not as exponential decay