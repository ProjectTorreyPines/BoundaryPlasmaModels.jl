function setup_model(boundary_plasma_model::LengyelModel, dd::IMAS.dd)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    hfs_sol2, lfs_sol2 = IMAS.sol(eqt, dd.wall, levels=2)
    sol = lfs_sol2[1]

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

    boundary_plasma_model.parameters.target.f_omp2target_expansion = (lfs_sol2[2].r[lfs_sol2[2].midplane_index] - lfs_sol2[1].r[lfs_sol2[1].midplane_index]) / sqrt((lfs_sol2[2].r[end] - lfs_sol2[1].r[end])^2 + (lfs_sol2[2].z[end] - lfs_sol2[1].z[end])^2)
    boundary_plasma_model.parameters.target.λ_target = boundary_plasma_model.parameters.sol.λ_omp * boundary_plasma_model.parameters.target.f_omp2target_expansion
    boundary_plasma_model.parameters.target.α_sp = atan(sol.Bp[end] / sol.Bt[end])
    boundary_plasma_model.parameters.target.f_pol_projection = tan(boundary_plasma_model.parameters.target.α_sp) # can be calculated internally
    boundary_plasma_model.parameters.target.θ_sp = sol.strike_angles[end]
    boundary_plasma_model.parameters.target.f_perp_projection = 1.0 / cos(boundary_plasma_model.parameters.target.θ_sp) # can be calculated internally
    boundary_plasma_model.parameters.target.f_spread_pfr = 1.0

    boundary_plasma_model.parameters.integral.T_down = 0.0
    boundary_plasma_model.parameters.integral.Zeff_exp = -0.3
    boundary_plasma_model.parameters.integral.Texp = 0.5
    boundary_plasma_model.parameters.integral.Lexp = 1.0
    boundary_plasma_model.parameters.integral.κ0 = 2390.0
end
