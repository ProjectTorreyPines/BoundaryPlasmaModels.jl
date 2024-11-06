#=
Author: Jerome Guterl (guterlj@fusion.gat.com)
 Company: General Atomics
 federoczak_peret_lambda_q.jl (c) 2024=#

mu_0 = 4e-7 * pi
mp = 1.67e-27
qe  = 1.6e-19

function compute_ρₛ(Te,A,Z,B)
    return sqrt(2*Te* A*mp/qe)/(B*Z)
end

function compute_σ_N(l_para, rhos)
    return 2*rhos/l_para
end

function compute_σ_ϕ(l_para, rhos)
    return 2*rhos/l_para
end

function compute_g(rhos,R,G0)
    return rhos/R*G0  
end

function compute_λ₀(g, sig_N, sig_phi, rhos, gamma)
    # alpha_ln = 3.9 
    # ln_names = ['g', 'sig_N', 'sig_PHI', 'rhos']
    # ln_powers = [3/11, -4/11, -2/11, 1]
    return 3.9 * g^(3/11) *sig_N^(-4/11) *sig_phi^(-2/11) * rhos * gamma^(-4/11)
end
function compute_λq(B_omp, Tesep, R_omp, l_para; G=1, α=3, γ=5.3,A=2,Z=1)
    ρₛ = compute_ρₛ(Tesep,A,Z,B_omp)
    g = compute_g(ρₛ,R_omp,G)
    σ_N = compute_σ_N(l_para, ρₛ)
    σ_ϕ = compute_σ_ϕ(l_para, ρₛ)
    λ₀ = compute_λ₀(g, σ_N,  σ_ϕ, ρₛ, γ)
    λₚ = compute_λₚ(λ₀, α)
    λq = compute_λq(λₚ,γ)
    return λq
end


function compute_λq(λₚ,γ)
    return 2/7*λₚ * sqrt(γ) / (sqrt(γ)-1)
end 

function compute_λₚ(λ₀, α)
    return λ₀ / α
end