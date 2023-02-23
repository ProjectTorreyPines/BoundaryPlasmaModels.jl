compute_P_SOL(dd::IMAS.dd) = compute_P_SOL( dd.core_profiles.profiles_1d[], dd.core_sources)
compute_P_SOL(cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources) = IMAS.total_sources(core_sources, cp1d).electrons.power_inside[end] + IMAS.total_sources(core_sources, cp1d).total_ion_power_inside[end]

compute_ϵ(dd) = compute_ϵ(dd.equilibrium.time_slice[])
compute_ϵ(eqt) = eqt.boundary.minor_radius / eqt.boundary.geometric_axis.r

get_a(dd) = get_a(dd.equilibrium.time_slice[])
get_a(eqt) = eqt.boundary.minor_radius

compute_mean_Bpol(dd) = compute_mean_Bpol(dd.equilibrium.time_slice[])
compute_mean_Bpol(eqt) = compute_mean_Bpol(get_a(eqt),get_κ(eqt),get_Ip(eqt))
compute_mean_Bpol(a::T, κ::T, Ip::T) where {T<:Real} = (constants.μ_0 * Ip) / (2π * a * sqrt((1.0 + κ^2) / 2.0))

get_κ(dd) = get_κ(dd.equilibrium.time_slice[])
get_κ(eqt) = eqt.profiles_1d.elongation[end]

get_Ip(dd) = get_Ip(dd.equilibrium.time_slice[])
get_Ip(eqt) = eqt.global_quantities.ip

get_Ip(dd) = get_Ip(dd.equilibrium.time_slice[])
get_Ip(eqt) = eqt.global_quantities.ip

compute_mean_R0(dd) = compute_mean_R0(dd.equilibrium.time_slice[])
compute_mean_R0(eqt) = (eqt.profiles_1d.r_outboard[end] .+ eqt.profiles_1d.r_inboard[end]) / 2.0

get_R0(dd) =  get_R0(dd.equilibrium.time_slice[])
get_R0(eqt) = eqt.global_quantities.magnetic_axis.r

