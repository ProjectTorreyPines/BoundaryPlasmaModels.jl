"test 
$(TYPEDSIGNATURES)
"
function test(aaaa::Float64,bbbb::Int64; cccc::Symbol = :testing)
    println(aaaa)
end
# using PlasmaFacingSurfaces
# import Term.Trees: Tree
# import PlasmaFacingSurfaces: get_projected_λ_imp, get_projected_λ_omp, Eqt,get_Bp_omp, PFSDesign, DivertorLeg, get_R_omp, get_R_imp,get_α_imp,get_α_omp,get_dps,get_α
# """Eich scaling (NF 53 093031)"""
# function λq_eich(R::T,P_SOL::T,Bpol::T,ϵ::T) where {T<:Float64}
  
#   if Bpol >0.25
#     error("Bpol>0.25T: no need to use the aspect ratio parameter... use λq_eich(Bpol::T)")
#     end
#     return 1.35e-3 * R^0.04 * Bpol^(-0.92) * ϵ^(0.42)
# end
# λq_eich(eqt::Eqt) = λq_eich(get_Bp_omp(eqt)) 

# λq_eich(Bpol::Float64) = 0.63e-3 * Bpol^(-1.19)


# # information
# function summary(P_SOL::Float64, pfs::PFSDesign; method_λ_q_omp = λq_eich, method_λ_q_imp = λq_eich, upper_lower_split::Float64 = 0.8, outer_inner_split::Float64 = 0.8)
#     eqt = pfs.eqt
#     @assert 0.0 <= upper_lower_split <= 1.0 
#     @assert 0.0 <= outer_inner_split <= 1.0 
#     info = Dict{Symbol,Any}() 
#     info[:P_SOL] = P_SOL
#     info[:outer_inner_split] = outer_inner_split
#     info[:upper_lower_split] = upper_lower_split
#     R_omp = get_R_omp(eqt)
#     R_imp = get_R_imp(eqt)
#     α_omp = get_α_omp(eqt)
#     α_imp = get_α_imp(eqt)
#     λ_q_omp = method_λ_q_omp(eqt)
#     λ_q_imp = method_λ_q_imp(eqt)
#     P_SOL_outer = P_SOL * outer_inner_split
#     P_SOL_inner = P_SOL * (1.0 -outer_inner_split)
#     Q_perp_omp = P_SOL_outer /(2π*R_omp*λ_q_omp)
#     Q_perp_imp = P_SOL_inner/(2π*R_imp*λ_q_imp)

#     Q_para_omp =  Q_perp_omp /sin(α_omp)
#     Q_para_imp =  Q_perp_imp /sin(α_imp)

# info[:outer] = Dict{Symbol,Any}()
# info[:inner] = Dict{Symbol,Any}()
# info[:outer][:P_in] = P_SOL_outer 
# info[:inner][:P_in] =  P_SOL_inner
# info[:P_SOL] = P_SOL
# info[:outer][:λ_q_mp] = λ_q_omp
# info[:inner][:λ_q_mp] = λ_q_imp
# info[:outer][:α_mp] = α_omp * 180/pi
# info[:inner][:α_mp] = α_imp * 180/pi
# info[:outer][:Q_para_mp] = Q_para_omp
# info[:outer][:Q_perp_mp] = Q_perp_omp
# info[:inner][:Q_para_mp] = Q_para_imp
# info[:inner][:Q_perp_mp] = Q_perp_imp
# info[:outer][:R] = R_omp
# info[:inner][:R] = R_imp
# if pfs.divertors.upper !== nothing
#     set_target_info!(info[:inner], pfs.divertors.upper.inner, eqt)
#     set_target_info!(info[:outer], pfs.divertors.upper.outer, eqt)
# end
# if pfs.divertors.lower !== nothing    
#     set_target_info!(info[:inner], pfs.divertors.lower.inner, eqt)
#     set_target_info!(info[:outer], pfs.divertors.lower.outer, eqt)
# end

# return info
# end 

# function set_target_info!(info_FS::Dict,leg::DivertorLeg, eqt::Eqt)
    
#     leg_id = Symbol(join(String.(leg.id), "_")) 
#     info_FS[leg_id] = Dict{Symbol,Any}()
#     info_FS[leg_id][:α] = get_α(get_dps(leg.target,:sp), eqt)[1] * 180/pi
#     info_FS[leg_id][:θ] = leg.target.θ * 180/pi
#     if :outer ∈ leg.id
#         info_FS[leg_id][:λ_projected] = get_projected_λ_omp(leg, info_FS[:λ_q_mp], eqt)
#     else
#         info_FS[leg_id][:λ_projected] = get_projected_λ_imp(leg, info_FS[:λ_q_mp], eqt)
#     end
#     info_FS[leg_id][:f_exp] = info_FS[leg_id][:λ_projected]/info_FS[:λ_q_mp] 
#     info_FS[leg_id][:Q_para_target] = info_FS[leg_id][:f_exp] * info_FS[:Q_para_mp]
#     info_FS[leg_id][:Q_perp_target] = info_FS[leg_id][:Q_para_target] * sind(info_FS[leg_id][:α])  
# end




# Eich scaling (NF 53 093031) without PSOL (because it is weak)
