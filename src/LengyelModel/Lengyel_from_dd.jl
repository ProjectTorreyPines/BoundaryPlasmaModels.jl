function setup_model(
    boundary_plasma_model::LengyelModel,
    dd::IMAS.dd;
    imp::Vector{Symbol}=[:Ne],
    f_imp::Vector{<:Real}=[0.02],
    f_spread_pfr::Real=1.0)

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    hfs_sol, lfs_sol = IMAS.sol(eqt, dd.wall, levels=2)
    P_SOL = IMAS.power_sol(dd.core_sources, cp1d)
    λ_omp = IMAS.widthSOL_eich(eqt, cp1d, dd.core_sources)

    setup_model(boundary_plasma_model, eqt, cp1d, lfs_sol, P_SOL, λ_omp; strike_index=0, imp, f_imp, f_spread_pfr)
end

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
* 0: automatic strike point selection favoring the outermost strike point
* 1: strike point at the 1st location in the `sol` OpenFieldLine
* anything else: strike point at the `end` location in the `sol` OpenFieldLine
"""
function setup_model(
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

    sol1 = sol[1]
    sol2 = sol[2]

    boundary_plasma_model.parameters.plasma.P_SOL = P_SOL
    boundary_plasma_model.parameters.plasma.R_omp = sol1.r[sol1.midplane_index]
    boundary_plasma_model.parameters.plasma.Ip = eqt.global_quantities.ip
    boundary_plasma_model.parameters.plasma.κ = eqt.boundary.elongation
    boundary_plasma_model.parameters.plasma.ϵ = eqt.boundary.minor_radius / eqt.boundary.geometric_axis.r
    boundary_plasma_model.parameters.plasma.Bt_omp = sol1.Bt[sol1.midplane_index]
    boundary_plasma_model.parameters.plasma.Bpol_omp = sol1.Bp[sol1.midplane_index]

    boundary_plasma_model.parameters.sol.n_up = cp1d.electrons.density_thermal[end]
    boundary_plasma_model.parameters.sol.T_up = cp1d.electrons.temperature[end]
    boundary_plasma_model.parameters.sol.λ_omp = λ_omp
    boundary_plasma_model.parameters.sol.imp = imp
    boundary_plasma_model.parameters.sol.f_imp = f_imp

    # automatic strike point selection (strike_index = 0) favors the outermost strike point
    if (strike_index == 1) || ((strike_index == 0) && (sol2.r[1] > sol2.r[end]))
        strike_index1 = 1
        strike_index2 = 1
    else
        strike_index1 = length(sol1.r)
        strike_index2 = length(sol2.r)
    end

    boundary_plasma_model.parameters.target.f_omp2target_expansion = sqrt((sol2.r[strike_index2] - sol1.r[strike_index1])^2 + (sol2.z[strike_index2] - sol1.z[strike_index1])^2) / (sol2.r[sol2.midplane_index] - sol1.r[sol1.midplane_index])
    boundary_plasma_model.parameters.target.f_spread_pfr = f_spread_pfr
    boundary_plasma_model.parameters.target.α_sp = atan(sol1.Bp[strike_index1] / sol1.Bt[strike_index1])
    boundary_plasma_model.parameters.target.θ_sp = sol1.strike_angles[strike_index1 == 1 ? 1 : 2] # [CURRENTLY NOT USED]
end

# questions:
# NO R_target ?
# Heat flux treated not as exponential decay