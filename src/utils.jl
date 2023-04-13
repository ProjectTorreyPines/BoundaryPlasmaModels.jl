compute_P_SOL(dd::IMAS.dd) = compute_P_SOL( dd.core_profiles.profiles_1d[], dd.core_sources)
compute_P_SOL(cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources) = IMAS.total_sources(core_sources, cp1d).electrons.power_inside[end] + IMAS.total_sources(core_sources, cp1d).total_ion_power_inside[end]

compute_ϵ(dd::IMAS.dd) = compute_ϵ(dd.equilibrium.time_slice[])
compute_ϵ(eqt::IMAS.equilibrium__time_slice) = eqt.boundary.minor_radius / eqt.boundary.geometric_axis.r

get_a(dd::IMAS.dd) = get_a(dd.equilibrium.time_slice[])
get_a(eqt::IMAS.equilibrium__time_slice) = eqt.boundary.minor_radius

compute_mean_Bpol(dd::IMAS.dd) = compute_mean_Bpol(dd.equilibrium.time_slice[])
compute_mean_Bpol(eqt::IMAS.equilibrium__time_slice) = compute_mean_Bpol(get_a(eqt),get_κ(eqt),get_Ip(eqt))
compute_mean_Bpol(a::T, κ::T, Ip::T) where {T<:Real} = (constants.μ_0 * Ip) / (2π * a * sqrt((1.0 + κ^2) / 2.0))

get_κ(dd::IMAS.dd) = get_κ(dd.equilibrium.time_slice[])
get_κ(eqt::IMAS.equilibrium__time_slice) = eqt.profiles_1d.elongation[end]

get_ϵ(dd::IMAS.dd) = get_ϵ(dd.equilibrium.time_slice[])
get_ϵ(eqt::IMAS.equilibrium__time_slice) = eqt.boundary.minor_radius/eqt.boundary.geometric_axis.r


get_Ip(dd::IMAS.dd) = get_Ip(dd.equilibrium.time_slice[])
get_Ip(eqt::IMAS.equilibrium__time_slice) = eqt.global_quantities.ip

compute_mean_R0(dd::IMAS.dd) = compute_mean_R0(dd.equilibrium.time_slice[])
compute_mean_R0(eqt::IMAS.equilibrium__time_slice) = (eqt.profiles_1d.r_outboard[end] .+ eqt.profiles_1d.r_inboard[end]) / 2.0

get_R0(dd::IMAS.dd) =  get_R0(dd.equilibrium.time_slice[])
get_R0(eqt::IMAS.equilibrium__time_slice) = eqt.global_quantities.magnetic_axis.r

get_PSOL(dd::IMAS.dd) = get_PSOL(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[], dd.core_sources)





function get_PSOL(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
    tot_src = IMAS.total_sources(core_sources, cp1d)
    return tot_src.electrons.power_inside[end] + tot_src.total_ion_power_inside[end]
end

get_ndiv(dd) = 1

get_omp_density(dd::IMAS.dd) = IMAS.interp1d(dd.core_profiles.profiles_1d[].grid.rho_tor_norm, dd.core_profiles.profiles_1d[].electrons.density_thermal).(1.00)
get_omp_Te(dd::IMAS.dd) = IMAS.interp1d(dd.core_profiles.profiles_1d[].grid.rho_tor_norm, dd.core_profiles.profiles_1d[].electrons.temperature).(1.00)

function get_target_outline(dd::IMAS.dd, location)
    target = get_target(dd,location)
    @assert target !== nothing "cannot find a target at location: $location "
    return FusionGeometryTools.RZLine(target.outline.r,target.outline.z) 
end

get_id(t::Union{IMAS.divertors__divertor___target,IMAS.equilibrium__time_slice___boundary_separatrix__strike_point}) = Symbol.(split(t.name))

function get_target(dd::IMAS.dd, location)
    for div in dd.divertors.divertor
        for t in div.target
            if all([l ∈ get_id(t) for l in location])
                return t
            end
        end
    end
    nothing
end
FusionGeometryTools.RZPoint(sp::IMASDD.equilibrium__time_slice___boundary_separatrix__strike_point) = FusionGeometryTools.RZPoint(sp.r,sp.z,get_id(sp))
function get_sp(dd::IMAS.dd, location::Vector{Symbol})
    sp = []
        for t in dd.equilibrium.time_slice[].boundary_separatrix.strike_point
            if all([l ∈ get_id(t) for l in location])
                push!(sp,FusionGeometryTools.RZPoint(t))
            end
        end
        if length(sp) ==  0 
            error("could not find a strike point in dd with location: $location")
        elseif length(sp) > 1
            error("multiple strike points found in dd with location: $location \n sp: $sp ")
        else
            return sp[1]
        end
end

function get_α_sp(dd :: IMAS.dd, location::Vector{Symbol}) 
    sp = get_sp(dd, location)
    return FusionGeometryTools.get_α(sp,dd.equilibrium) 
end

function get_θ_sp(dd ::IMAS.dd, location::Vector{Symbol}) 
    target = get_target(dd::IMAS.dd, location)
    return target.angle_sp
end

function get_projected_λ_omp_target(λ_omp::Float64,dd::IMAS.dd; location=[:outer,:upper])
    target_outline = get_target_outline(dd, location)
    sp = get_sp(dd,location)
    p = FusionGeometryTools.get_omp_projection_target(λ_omp,target_outline,dd.equilibrium.time_slice[])
    @assert sp !==nothing
    FusionGeometryTools.curvilinear_distance(target_outline,sp,p)
end

function get_projected_flux_expension(λ_omp::Float64,dd::IMAS.dd; location=[:outer,:upper])
    return get_projected_λ_omp_target(λ_omp,dd; location)/λ_omp
end



