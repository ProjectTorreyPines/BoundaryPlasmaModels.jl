function setup_model(boundary_plasma_model::LengyelModel, dd::IMAS.dd)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    hfs_sol, lfs_sol = IMAS.sol(eqt, dd.wall, levels=2)
    sol = lfs_sol[1]
    sol2 = lfs_sol[2]

    boundary_plasma_model.parameters.plasma.P_SOL = IMAS.power_sol(dd.core_sources, cp1d)
    boundary_plasma_model.parameters.plasma.R_omp = sol.r[sol.midplane_index]
    boundary_plasma_model.parameters.plasma.Ip = eqt.global_quantities.ip
    boundary_plasma_model.parameters.plasma.κ = eqt.boundary.elongation
    boundary_plasma_model.parameters.plasma.ϵ = eqt.boundary.minor_radius / eqt.boundary.geometric_axis.r
    boundary_plasma_model.parameters.plasma.Bt_omp = sol.Bt[sol.midplane_index]
    boundary_plasma_model.parameters.plasma.Bpol_omp = sol.Bp[sol.midplane_index]

    boundary_plasma_model.parameters.sol.n_up = cp1d.electrons.density_thermal[end]
    boundary_plasma_model.parameters.sol.T_up = cp1d.electrons.temperature[end]

    boundary_plasma_model.parameters.sol.λ_omp = IMAS.widthSOL_eich(eqt, cp1d, dd.core_sources)
    boundary_plasma_model.parameters.sol.f_imp = [0.02]
    boundary_plasma_model.parameters.sol.imp = [:Ne]

    boundary_plasma_model.parameters.target.f_omp2target_expansion = sqrt((sol2.r[end] - sol.r[end])^2 + (sol2.z[end] - sol.z[end])^2) / (sol2.r[sol2.midplane_index] - sol.r[sol.midplane_index])
    boundary_plasma_model.parameters.target.f_spread_pfr = 1.0
    boundary_plasma_model.parameters.target.α_sp = atan(sol.Bp[end] / sol.Bt[end])
    boundary_plasma_model.parameters.target.θ_sp = sol.strike_angles[end] # [CURRENTLY NOT USED]
end

# questions:
# NO R_target ?
# Heat flux treated not as exponential decay