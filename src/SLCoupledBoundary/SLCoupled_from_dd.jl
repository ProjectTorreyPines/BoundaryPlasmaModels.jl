"""
    setup_model(
        model::SLCoupledModel,                            
        λ_omp::Real,
        eqt::IMAS.equilibrium__time_slice,
        cp1d::IMAS.core_profiles__profiles_1d,
        cs::IMAS.core_sources,
        sol1::IMAS.OpenFieldLine,
        impurities::Vector{Symbol},
        impurities_fraction::Vector{<:Real},
        heat_spread_factor::Real = 1.0,
        strike_index::Int,
    )

Set up parameters for the SLCoupledModel using target-side inputs.

`strike_index` is the index of the strike point in the points of the `sol1` OpenFieldLine:
`1` (its first point) or `length(sol1.r)` (its last point)

Upstream (Te_up, ne_up) is solved later; here we only populate inputs/storage.
"""
function setup_model(
    model::SLCoupledModel,                            
    λ_omp::Real,
    eqt::IMAS.equilibrium__time_slice,
    cp1d::IMAS.core_profiles__profiles_1d,
    cs::IMAS.core_sources,
    sol1::IMAS.OpenFieldLine;
    impurities::Vector{Symbol},
    impurities_fraction::Vector{<:Real},
    heat_spread_factor::Real = 1.0,
    strike_index::Int,
)
    @assert length(impurities) == length(impurities_fraction)

    model.eqt    = eqt
    model.cp1d   = cp1d
    model.sol1   = sol1
    p = model.parameters

    # plasma 
    p.plasma.P_SOL     = (IMAS.power_sol(cs, cp1d) / 2.0) # two strike points
    p.plasma.R_omp     = sol1.r[sol1.midplane_index]
    p.plasma.R_x      = eqt.boundary.x_point[1].r
    p.plasma.Ip        = eqt.global_quantities.ip
    p.plasma.κ         = eqt.boundary.elongation
    p.plasma.ϵ         = eqt.boundary.minor_radius / eqt.boundary.geometric_axis.r
    p.plasma.Bt_omp    = sol1.Bt[sol1.midplane_index]
    p.plasma.Bpol_omp  = sol1.Bp[sol1.midplane_index]

    # sol
    p.sol.λ_omp = λ_omp
    p.sol.Te_up = cp1d.electrons.density_thermal[end]
    p.sol.ne_up = cp1d.electrons.temperature[end]

    # impurities for radiation model
    p.radiation.species   = collect(impurities)
    p.radiation.fractions = collect(float.(impurities_fraction))

    # target geometry, at the strike point this model is set up for
    @assert strike_index in (1, length(sol1.r)) "strike_index must be 1 or length(sol1.r)"
    p.target.f_omp2target_expansion = sol1.total_flux_expansion[strike_index]
    p.target.R_strike = sol1.r[strike_index]
    p.target.L_para   = sol1.s[strike_index]
    p.target.f_spread_pfr = heat_spread_factor
    p.target.α_sp = atan(sol1.Bp[strike_index] / sol1.Bt[strike_index])
    p.target.θ_sp = sol1.strike_angles[strike_index == 1 ? 1 : 2]

    return model
end

