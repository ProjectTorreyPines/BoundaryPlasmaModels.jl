module SLCoupledFullModel

import SimulationParameters              # for setup_parameters!(...)
import IMAS
using IMAS                               # IMAS types used by from_dd & handles
import ADAS                              # for Lengyel model
import ..BoundaryPlasmaModels: SOLBoundaryModel, SOLBoundaryParameters
using Format
using Optim

export SLCoupledParameters, SLCoupledModel, setup_model, summary
export solve_target_known!,solve_upstream_known!

include("SLCoupled_parameters.jl")

const eV_to_J = 1.602176634e-19         # [J/eV]
const m_p     = 1.67262192369e-27       # [kg] 

Base.@kwdef mutable struct SLCoupledResults{T}
    # geometry / projections
    λ_target::T              = zero(T)
    f_pol_projection::T      = zero(T)   # ≈ tan(α_sp)

    # heat fluxes
    q_parallel_target::T     = zero(T)   # q∥ at target (unprojected)
    q_parallel_upstream::T   = zero(T)   # q∥ at OMP/upstream
    q_perp_target::T         = zero(T)   # projected to target surface
    q_perp_target_spread::T  = zero(T)   # after spread factor
    q_rad::T                 = zero(T)   # integrated radiated heat flux density (input/est.)    
    f_cool::T                = zero(T)   # fitted cooling factor
    f_mom::T                 = zero(T)   # fitted momentum factor


    # thermodynamics
    Te_up::T                 = zero(T)
    ne_up::T                 = zero(T)
    Te_t::T                  = zero(T)
    ne_t::T                  = zero(T)

    # impurities/radiation diagnostics
    #Zeff::T                  = zero(T)
    cz_total::T              = zero(T)
    cz_per_species::Vector{T}= T[]
end

function Base.show(io::IO, r::SLCoupledResults)
    for f in propertynames(r)
        println(io, "$(f): ", getfield(r, f))
    end
end

# --- model type ---
mutable struct SLCoupledModel <: SOLBoundaryModel
    parameters::SLCoupledParameters
    results::SLCoupledResults
    target::Union{Nothing, IMAS.divertors__divertor___target}
    eqt::Union{Nothing, IMAS.equilibrium__time_slice}
    cp1d::Union{Nothing, IMAS.core_profiles__profiles_1d}
    sol1::Union{Nothing, IMAS.OpenFieldLine}
end

SLCoupledModel() = SLCoupledModel(
    SLCoupledParameters{Float64}(),
    SLCoupledResults{Float64}(),
    nothing, nothing, nothing, nothing
)
include("SLCoupled_from_dd.jl")
# ---------------------------
# helpers: polynomials & sheath
# ---------------------------

@inline function f_from_poly_Te(Te_t::Number, coeffs::NTuple{5,<:Real})
    ξ = log10(Te_t)
    y = coeffs[1] + coeffs[2]*ξ + coeffs[3]*ξ^2 + coeffs[4]*ξ^3 + coeffs[5]*ξ^4
    T10 = oftype(Te_t, 10.0)                  #
    f = one(Te_t) - (T10^y)                   
    return clamp(f, zero(Te_t), one(Te_t))
end

@inline γ_sh(τt::Number) = 5.69 + 3.0*τt + 0.5*log1p(τt)

@inline function c_s_t(Te_t::Number, τt::Number, mi_amu::Number)
    mi = mi_amu * m_p                        
    return sqrt(eV_to_J * Te_t * (one(Te_t)+τt) / mi)
end

@inline function q_parallel_at_target(Te_t::Number, ne_t::Number, τt::Number, mi_amu::Number)
    cs = c_s_t(Te_t, τt, mi_amu)
    return γ_sh(τt) * ne_t * (eV_to_J*Te_t) * cs
end

@inline sin_pitch(Bp::Number, Bt::Number) = Bp / sqrt(Bp^2 + Bt^2)


# ---------------------------
# ADAS LINT (experience-corrected, single mixture integral)
# LINT = ∫ [ T^Texp * Zeff_mix(fracs, ne, T)^Zeff_exp * (Σ_i f_i * Lz_i(ne, T))^Lexp ] dT
# default Texp=0.5 (≈√T), Lexp=1.0 (linear), Zeff_exp=-0.3 (empirical)
# ---------------------------
function _LINT_adas_sum(Tmin::Float64, Tmax::Float64,
                        species::Vector{Symbol}, fracs::Vector{Float64};
                        ne::Float64 = 1e20, N::Int = 400,
                        Zeff_exp::Float64 = -0.3,
                        Texp::Float64 = 0.5,
                        Lexp::Float64 = 1.0)

    # empty mix or degenerate interval → 0
    if isempty(species) || isempty(fracs) || Tmin == Tmax
        return 0.0
    end
    @assert length(species) == length(fracs) "species and fracs must have same length"

    # build temperature grid
    Tlo = min(Tmin, Tmax); Thi = max(Tmin, Tmax)
    Tgrid = collect(range(Tlo, Thi; length=N))
    Δ = (Thi - Tlo) / (N - 1)

    # mixture Zeff callable (handles all species together)
    zeff_mix = ADAS.get_Zeff(species)   # usage: zeff_mix(fracs, ne, T)

    # fetch cooling interpolants per species (prefer Lz_tot, fallback Lztot)
    Lz_funcs = Vector{Function}(undef, length(species))
    for (k, sym) in pairs(species)
        data = ADAS.get_cooling_rates(sym)
        Lz_funcs[k] = (n, Te) -> data.Lztot(n, Te)
    end

    # integrand values at each T: T^Texp * (Σ f_i Lz_i)^Lexp
    vals = Vector{Float64}(undef, N)
    @inbounds for (j, Te) in enumerate(Tgrid)
        # mixture Lz at (ne, Te)
        Lmix = 0.0
        @inbounds for (fi, Lz) in zip(fracs, Lz_funcs)
            fi == 0 && continue
            Lmix += fi * Lz(ne, Te)
        end
        # Zeff correction and temperature weighting
        #Zcorr = zeff_mix(fracs, ne, Te)
        vals[j] = (Te^Texp) * (Lmix^Lexp)
    end

    # trapezoidal integration over T
    return Δ * (0.5 * vals[1] + sum(@view vals[2:end-1]) + 0.5 * vals[end])
end

# =============== Residuals (R1, R2, R3) for optimization ===============
"""
Residual system using your three equations:

R1: (sheath heat flux * Rt) - ((1 - fcool(Te_t)) * q_leg * Ru) = 0
R2: nt*Tt*(1+τt)*(1+Mt^2) - nu*Tu*(1+τu)*(1+Mu^2)*(1 - fmom(Te_t)) = 0
R3: Tu^(7/2) - Tt^(7/2) - (7/2)*(Rt/Ru)*[ q_leg * f_cond * L_par ] / [ k_e * (1 - fcool(Te_t)) ] = 0

Here q_leg = P_SOL / (2π Ru λ_q sinθ), with sinθ = Bp / sqrt(Bp^2 + Bt^2).
"""
function residuals_vec(x; mode::Symbol, m::SLCoupledModel, do_print::Bool=false)
    x = abs.(x)
    x[3] = IMAS.mirror_bound(x[3], 0.0, 1.0)

    p = m.parameters

    τt, τu = p.coupling.tau_t, p.coupling.tau_u
    mi     = p.coupling.mi_amu
    Ru, Rt = p.plasma.R_omp, p.target.R_strike
    λq     = p.sol.λ_omp
    Bp, Bt = p.plasma.Bpol_omp, p.plasma.Bt_omp
    Lpar   = p.target.L_para
    ke     = p.coupling.kappa_e
    fcond  = p.coupling.f_cond
    Mt, Mu = p.coupling.M_t, p.coupling.M_u

    sθ   = sin_pitch(Bp, Bt)
    qleg = p.plasma.P_SOL / (2π * Ru * λq * sθ)   # [W/m^2] at OMP

    if mode === :target_known
        # unknowns: [Tu, nu, fc]; known: Te_t, ne_t
        Tu, nu, fc = x
        Tt, nt = p.target.Te_t, p.target.ne_t
    elseif mode === :upstream_known
        # unknowns: [Tt, nt, fc]; known: Te_u, ne_u
        Tt, nt, fc = x
        Tu, nu = p.sol.Te_up, p.sol.ne_up
    else
        error("mode must be :target_known or :upstream_known")
    end

    # Fractions derived from Te_t via stored polynomial coefficients
    fm = f_from_poly_Te(Tt, p.coupling.fmom_poly)

    # R1: sheath vs available power mapping (units → [W/m] on both sides)
    qt_sheath = q_parallel_at_target(Tt, nt, τt, mi)  # [W/m^2]
    R1 = qt_sheath * Rt - (1 - fc) * qleg * Ru
    R1 = R1 / qleg
  
    # R2: momentum/pressure balance with fmom(Te_t)
    R2 = nt*Tt*(1+τt)*(1+Mt^2) - nu*Tu*(1+τu)*(1+Mu^2) * (1 - fm)
    R2 = R2 / nt*Tt

    # R3: Spitzer–Härm conduction integral with geometry and f_cond
    R3 = Tu^(3.5) - Tt^(3.5) - (7/2) * (Rt/Ru) * (qleg * fcond * Lpar) / (ke * (1 - fc))
    R3 = R3 / Tt^(3.5)

    return (R1, R2, R3)
end

# =============== Post-processing (populate results similar to your original solver) ===============
function postprocess!(m::SLCoupledModel)
    p = m.parameters
    r = m.results

    Rt, Ru   = p.target.R_strike, p.plasma.R_omp
    λq       = p.sol.λ_omp
    Lpar     = p.target.L_para
    ke       = p.coupling.kappa_e
    fcond    = p.coupling.f_cond
    τt       = p.coupling.tau_t
    mi       = p.coupling.mi_amu
    Bp, Bt   = p.plasma.Bpol_omp, p.plasma.Bt_omp
    sθ       = sin_pitch(Bp, Bt)
    qleg     = p.plasma.P_SOL / (2π * Ru * λq * sθ)

    # Fractions at current Te_t
    #fc = f_from_poly_Te(p.target.Te_t, p.coupling.fcool_poly)
    fc = p.coupling.fcool_optim
    fm = f_from_poly_Te(p.target.Te_t, p.coupling.fmom_poly)
    r.f_cool = fc
    r.f_mom  = fm

    # Heat fluxes
    q_par_u = qleg 
    q_par_t = q_par_u * (1.0 - fc)
    r.q_parallel_target   = q_par_t
    r.q_parallel_upstream = q_par_u

    r.f_pol_projection    = tan(p.target.α_sp)
    r.λ_target            = p.sol.λ_omp * p.target.f_omp2target_expansion
    r.q_perp_target       = r.q_parallel_target * r.f_pol_projection
    r.q_perp_target_spread= r.q_perp_target / p.target.f_spread_pfr

    # Simple q_rad estimate (geometric mean form used previously)
    r.q_rad = sqrt(max(0.0, q_par_u^2 - q_par_t^2))

    # Radiation split to species (ADAS LINT integral)
    if !isempty(p.radiation.species) && !isempty(p.radiation.fractions)
        Tmin = Float64(min(p.target.Te_t, p.sol.Te_up))
        Tmax = Float64(max(p.target.Te_t, p.sol.Te_up))
        LINT = _LINT_adas_sum(
            Tmin, Tmax,
            Vector{Symbol}(p.radiation.species),
            Vector{Float64}(p.radiation.fractions);
            ne = Float64(p.sol.ne_up)
        )
        if LINT > 0
            cz_total = ((Float64(r.q_rad)^2) /
                (2 * Float64(ke) * Float64(p.sol.ne_up)^2 * Float64(p.sol.Te_up)^2 * LINT)) / 5.0
            w = copy(Vector{Float64}(p.radiation.fractions))
            s = max(sum(w), eps(Float64))
            w ./= s
            r.cz_per_species = cz_total .* w
            r.cz_total = sum(r.cz_per_species)
        else
            r.cz_per_species = Float64[]
            r.cz_total = 0.0
        end
    else
        r.cz_per_species = Float64[]
        r.cz_total = 0.0
    end

    # Mirror key state to results
    r.Te_up = p.sol.Te_up
    r.ne_up = p.sol.ne_up
    r.Te_t  = p.target.Te_t
    r.ne_t  = p.target.ne_t

    return m
end

# =============== Public solvers (optimization-only) ===============
function run!(m::SLCoupledModel; mode::Symbol=:target_known,
              x0=nothing, lb=nothing, ub=nothing, weights=(1.0,1.0,1.0),
              kwargs...)
    if mode === :target_known
        solve_target_known!(m; x0=x0, lb=lb, ub=ub, weights=weights)
    elseif mode === :upstream_known
        solve_upstream_known!(m; x0=x0, lb=lb, ub=ub, weights=weights)
    else
        error("Unsupported mode = $mode")
    end
    return m.results
end

function (m::SLCoupledModel)(; kwargs...)
    return run!(m; kwargs...)   
end
"""
    solve_target_known!(m; x0=nothing, lb=nothing, ub=nothing, weights=(1,1,1))

Target known: P_SOL, Te_t, ne_t must be set in `m.parameters`.
Solves for upstream (Te_up, ne_up) by minimizing the weighted sum of squared residuals.
Writes the solution into `m.parameters` and updates `m.results` via `postprocess!`.
Returns a named tuple with the optimizer output for diagnostics.
"""
function solve_target_known!(m::SLCoupledModel; x0=nothing, lb=nothing, ub=nothing, weights=(1.0,1.0,1.0))
    p = m.parameters

    cost(x; do_print=false) = begin
        #push!(tmp_in, x)
        r = residuals_vec([x[1],x[2] * 1E20, x[3]]; mode=:target_known, m, do_print)
        rw = r .* weights
        c = sqrt(sum(rw.^2)) 
        return c
    end

    # Unknowns x = [Tu, nu, fc]
    # density in units of 1E20
    x0 === nothing && (x0 = [p.target.Te_t*5, p.target.ne_t/1E20/5, 0.9])
    lb === nothing && (lb = [0.1, 1E-3, 1E-3])
    ub === nothing && (ub = [1e4, 10.0, 1.0 - 1E-3])

    #res = optimize(cost, lb, ub, x0, Fminbox(LBFGS()))

    # rough initial guess of x0
    factors = range(1.0, 20.0, 10)
    costs = [cost([p.target.Te_t*factor, p.target.ne_t/1E20/factor, 0.9]) for factor in factors]
    best_factor = factors[argmin(costs)]

    # run actual optimization
    x0 = [p.target.Te_t*best_factor, p.target.ne_t/1E20/best_factor, 0.9]
    res = optimize(cost, x0, Newton(); autodiff=:forward)

    Tu, nu, fc = abs.(res.minimizer)
    fc = IMAS.mirror_bound(fc, 0.0, 1.0)
    p.sol.Te_up = Tu
    p.sol.ne_up = nu * 1E20
    p.coupling.fcool_optim = fc

    postprocess!(m)
    Tu, nu_s, fc = abs.(res.minimizer)         
    println("Final minimizer: ", Float64.(abs.(res.minimizer)))
    println("Final residuals: ",Float64.(residuals_vec([Tu, nu_s*1e20, fc]; mode=:target_known, m=m)))
    println("Final cost: ", res.minimum)
    return res
end

# """
#     solve_upstream_known!(m; x0=nothing, lb=nothing, ub=nothing, weights=(1,1,1))

# Upstream known: P_SOL, Te_up, ne_up must be set in `m.parameters`.
# Solves for target (Te_t, ne_t) by minimizing the weighted sum of squared residuals.
# Writes the solution into `m.parameters` and updates `m.results` via `postprocess!`.
# Returns a named tuple with the optimizer output for diagnostics.
# """
function solve_upstream_known!(m::SLCoupledModel; x0=nothing, lb=nothing, ub=nothing, weights=(1.0,1.0,1.0))
    p = m.parameters

    cost(x; do_print=false) = begin
        r = residuals_vec([x[1], x[2] * 1E20, x[3]]; mode=:upstream_known, m, do_print)
        rw = r .* weights
        c = sqrt(sum(rw.^2)) 
        return c
    end

    # Unknowns x = [Tt, nt, fc]
    x0 === nothing && (x0 = [p.sol.Te_up/5, p.sol.ne_up/1E20*5, 0.9])
    lb === nothing && (lb = [0.1, 1E-3, 1E-3])
    ub === nothing && (ub = [1e3, 10.0, 1.0 - 1E-3])

    # rough initial guess of x0
    factors = range(1.0, 20.0, 10)
    costs = [cost([p.sol.Te_up*factor, p.sol.ne_up/1E20/factor, 0.9]) for factor in factors]
    best_factor = factors[argmin(costs)]

    x0 = [p.sol.Te_up*best_factor, p.sol.ne_up/1E20/best_factor, 0.9]
    res = optimize(cost, x0, Newton(); autodiff=:forward)

    Tt, nt, fc = abs.(res.minimizer)
    fc = IMAS.mirror_bound(fc, 0.0, 1.0)
    p.target.Te_t = Tt
    p.target.ne_t = nt * 1E20
    p.coupling.fcool_optim = fc

    postprocess!(m)

    Tt, nt_s, fc = abs.(res.minimizer)
    println("Final minimizer: ", Float64.(abs.(res.minimizer)))
    println("Final residuals: ", Float64.(residuals_vec([Tt, nt_s*1e20, fc]; mode=:upstream_known, m=m)))
    println("Final cost: ", res.minimum)

    return res
end

# =============== Summary printer ===============
function Base.summary(m::SLCoupledModel)
    p, r = m.parameters, m.results

    printfmtln("SLCoupled (Optimization-based, 3-equation system)")
    println("Upstream (solved or given)")
    printfmtln("├─ {:<22} = {:.2f} eV",  "Te_up", r.Te_up)
    printfmtln("├─ {:<22} = {:.3e} m^-3","ne_up", r.ne_up)
    printfmtln("├─ {:<22} = {:.3f} m",   "λ_omp", p.sol.λ_omp)
    printfmtln("├─ {:<22} = {:.3f} m",   "R_omp", p.plasma.R_omp)
    printfmtln("├─ {:<22} = {:.2f} T",   "Bp_omp", p.plasma.Bpol_omp)
    printfmtln("└─ {:<22} = {:.2f} T",   "Bt_omp", p.plasma.Bt_omp)
    println()

    println("Target (solved or given)")
    printfmtln("├─ {:<22} = {:.2f} eV",  "Te_t", r.Te_t)
    printfmtln("├─ {:<22} = {:.3e} m^-3","ne_t", r.ne_t)
    printfmtln("├─ {:<22} = {:.3f} m",   "L_para", p.target.L_para)
    printfmtln("├─ {:<22} = {:.3f}",     "f_omp→target", p.target.f_omp2target_expansion)
    printfmtln("├─ {:<22} = {:.3f} m",   "R_strike", p.target.R_strike)
    printfmtln("├─ {:<22} = {:.1f} deg", "α_pitch", p.target.α_sp*180/π)
    printfmtln("└─ {:<22} = {:.1f} deg", "θ_target", p.target.θ_sp*180/π)
    println()

    println("Coupling (evaluated at Te_t)")
    printfmtln("├─ {:<22} = {:.3f}", "f_cool(Te_t)", r.fcool_optim)
    printfmtln("└─ {:<22} = {:.3f}", "f_mom(Te_t)",  r.f_mom)
    println()

    println("Radiation / heat flux")
    printfmtln("├─ {:<22} = {:.3e} W/m²", "q_parallel (target)", r.q_parallel_target)
    printfmtln("├─ {:<22} = {:.3e} W/m²", "q_parallel (upstream)", r.q_parallel_upstream)
    printfmtln("├─ {:<22} = {:.3e} W/m²", "q_perp target", r.q_perp_target)
    printfmtln("├─ {:<22} = {:.3e} W/m²", "q_perp target (spread)", r.q_perp_target_spread)
    printfmtln("└─ {:<22} = {:.3e} W/m²", "q_rad (est.)", r.q_rad)
    println()

    if !isempty(r.cz_per_species)
        pairstr = join(["$(p.radiation.species[i])=$(r.cz_per_species[i])" for i in eachindex(p.radiation.species)], ", ")
        printfmtln("cz per species: {}", pairstr)
        printfmtln("cz total     : {:.3e}", r.cz_total)
    else
        println("cz per species: (n/a)")
    end
    nothing
end

end # module