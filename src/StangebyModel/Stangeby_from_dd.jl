"""
    setup_model(
      boundary_plasma_model::StangebyModel,
      λ_omp::Real,
      eqt::IMAS.equilibrium__time_slice,
      cp1d::IMAS.core_profiles__profiles_1d,
      cs::IMAS.core_sources,
      sol1::IMAS.OpenFieldLine;
      impurities::Vector{Symbol},
      impurities_fraction::Vector{<:Real},
      heat_spread_factor::Real=1.0,
      strike_index::Int)

`strike_index` is the index of the strike point in the points of the `sol1` OpenFieldLine:
`1` (its first point) or `length(sol1.r)` (its last point)
"""
function setup_model(
  boundary_plasma_model::StangebyModel,
  λ_omp::Real,
  eqt::IMAS.equilibrium__time_slice,
  cp1d::IMAS.core_profiles__profiles_1d,
  cs::IMAS.core_sources,
  sol1::IMAS.OpenFieldLine;
  impurities::Vector{Symbol},
  impurities_fraction::Vector{<:Real},
  heat_spread_factor::Real=1.0,
  strike_index::Int)

    @assert (length(impurities) == length(impurities_fraction))
    @assert strike_index in (1, length(sol1.r)) "strike_index must be 1 or length(sol1.r)"
    n_targets = 2
    boundary_plasma_model.parameters.plasma.P_SOL = IMAS.power_sol(cs, cp1d) / n_targets
    boundary_plasma_model.parameters.plasma.R_omp = sol1.r[sol1.midplane_index]
    boundary_plasma_model.parameters.plasma.Ip = eqt.global_quantities.ip
    boundary_plasma_model.parameters.plasma.κ = eqt.boundary.elongation
    boundary_plasma_model.parameters.plasma.ϵ = eqt.boundary.minor_radius / eqt.boundary.geometric_axis.r
    boundary_plasma_model.parameters.plasma.Bt_omp = sol1.Bt[sol1.midplane_index]
    boundary_plasma_model.parameters.plasma.Bpol_omp = sol1.Bp[sol1.midplane_index]

    boundary_plasma_model.parameters.sol.T_up = cp1d.electrons.temperature[end]
    boundary_plasma_model.parameters.sol.λ_omp = λ_omp

    boundary_plasma_model.parameters.target.f_omp2target_expansion = sol1.total_flux_expansion[strike_index]
    λ_target = boundary_plasma_model.parameters.sol.λ_omp * boundary_plasma_model.parameters.target.f_omp2target_expansion
    boundary_plasma_model.parameters.target.R = sol1.r[strike_index]
    boundary_plasma_model.parameters.target.L_para = sol1.s[strike_index]
    boundary_plasma_model.parameters.target.f_spread_pfr = heat_spread_factor
    boundary_plasma_model.parameters.target.α_sp = atan(sol1.Bp[strike_index] / sol1.Bt[strike_index])
    return boundary_plasma_model.parameters.target.θ_sp = sol1.strike_angles[strike_index == 1 ? 1 : 2]
end