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
    # Cap Te_t at 10.0
    Te_eff = clamp(Te_t, oftype(Te_t, 0.1), oftype(Te_t, 10.0))
    ξ = log10(Te_eff)
    y = coeffs[1] + coeffs[2]*ξ + coeffs[3]*ξ^2 + coeffs[4]*ξ^3 + coeffs[5]*ξ^4
    T10 = oftype(Te_t, 10.0)                  #
    f = one(Te_t) - (T10^y)                  
    ϵ = oftype(Te_t, 1e-6)
    return clamp(f, ϵ, one(Te_t) - ϵ)
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

# =============== Post-processing (populate results) ===============
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
    Psol_eff = p.plasma.P_SOL * p.plasma.f_psol 
    qleg     = Psol_eff / (2π * Ru * λq * sθ)

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
function run!(m::SLCoupledModel; mode::Symbol=:target_known, x0=nothing, kwargs...)
    if mode === :target_known
        solve_target_known!(m; x0=x0)
    elseif mode === :upstream_known
        solve_upstream_known!(m; x0=x0)
    else
        error("Unsupported mode = $mode")
    end
    return m.results
end

function (m::SLCoupledModel)(; kwargs...)
    return run!(m; kwargs...)   
end
"""
    solve_target_known!(m; x0=nothing)

 Direct triangular solve (no optimizer):
- R1 → compute `f_cool` from sheath vs. available power
- R3 → compute upstream temperature `Te_up`
- R2 → compute upstream density `ne_up`

Inputs required in `m.parameters`: `P_SOL`, target `(Te_t, ne_t)`, geometry, and coupling constants.
Writes the solution into `m.parameters` and updates `m.results` via `postprocess!`.
Returns a small named tuple mimicking optimizer results for compatibility.
"""
function solve_target_known!(m::SLCoupledModel; x0=nothing)  
    p = m.parameters
    ϵ = 1e-6

    # read constants/parameters
    τt, τu = p.coupling.tau_t, p.coupling.tau_u
    mi     = p.coupling.mi_amu
    Ru, Rt = p.plasma.R_omp, p.target.R_strike
    λq     = p.sol.λ_omp
    Bp, Bt = p.plasma.Bpol_omp, p.plasma.Bt_omp
    Lpar   = p.target.L_para
    ke     = p.coupling.kappa_e
    fcond  = p.coupling.f_cond
    Mt, Mu = p.coupling.M_t, p.coupling.M_u

    # known target conditions
    Tt, nt = p.target.Te_t, p.target.ne_t

    # q_leg
    sθ   = max(sin_pitch(Bp, Bt), 1e-12)
    Psol_eff = p.plasma.P_SOL * p.plasma.f_psol 
    qleg = Psol_eff / (2π * Ru * λq * sθ)   # [W/m^2]

    # R1: f_cool
    qt  = q_parallel_at_target(Tt, nt, τt, mi)    # [W/m^2]
    fc  = 1 - (qt * Rt) / (qleg * Ru)
    fc  = clamp(fc, ϵ, 1 - ϵ)

    # R3: Te_up
    Tu = (Tt^(3.5) + (7/2) * (Rt/Ru) * (qleg * fcond * Lpar) / (ke * (1 - fc)))^(2/7)

    # R2: ne_up
    fm = f_from_poly_Te(Tt, p.coupling.fmom_poly)
    nu = nt * Tt * (1 + τt) * (1 + Mt^2) / ( Tu * (1 + τu) * (1 + Mu^2) * (1 - fm) )

    # write back & postprocess
    p.sol.Te_up = Tu
    p.sol.ne_up = nu
    p.coupling.fcool_optim = fc
    postprocess!(m)

    println("Direct solve (target_known) → Tu=$(Float64(Tu)), nu=$(Float64(nu)), fc=$(Float64(fc))") 

    # return an optimizer-like result for compatibility
    return (minimizer = [Tu, nu/1e20, fc], minimum = 0.0, converged = true)  
end

"""
    solve_upstream_known!(m; x0=nothing)

 Direct 1D solve (Brent) when upstream is known:
- Use R2 to express `n_t(T_t)`
- Use R1 to express `f_cool(T_t)`
- Insert into R3 to define a scalar residual `g(T_t)` and minimize `g(T_t)^2` on [T_min, 0.99*Te_up]

Inputs required in `m.parameters`: `P_SOL`, upstream `(Te_up, ne_up)`, geometry, and coupling constants.
Writes the solution into `m.parameters` and updates `m.results` via `postprocess!`.
Returns a small named tuple mimicking optimizer results for compatibility.
"""
function solve_upstream_known!(m::SLCoupledModel; x0=nothing)  
    p = m.parameters
    ϵ = 1e-6

    # constants/parameters
    Tu, nu = p.sol.Te_up, p.sol.ne_up
    τt, τu  = p.coupling.tau_t, p.coupling.tau_u
    mi      = p.coupling.mi_amu
    Ru, Rt  = p.plasma.R_omp, p.target.R_strike
    λq      = p.sol.λ_omp
    Bp, Bt  = p.plasma.Bpol_omp, p.plasma.Bt_omp
    Lpar    = p.target.L_para
    ke      = p.coupling.kappa_e
    fcond   = p.coupling.f_cond
    Mt, Mu  = p.coupling.M_t, p.coupling.M_u

    sθ   = max(sin_pitch(Bp, Bt), 1e-12)
    Psol_eff = p.plasma.P_SOL * p.plasma.f_psol 
    qleg = Psol_eff / (2π * Ru * λq * sθ)

    # R2 ⇒ n_t(T_t)
    C = nu*Tu*(1+τu)*(1+Mu^2) / ((1+τt)*(1+Mt^2))
    nt_of = Tt -> max( C * (1 - f_from_poly_Te(Tt, p.coupling.fmom_poly)) / max(Tt, ϵ), ϵ )

    # R1 ⇒ f_cool(T_t)
    qt_sheath = (Tt, nt) -> q_parallel_at_target(Tt, nt, τt, mi)
    fc_of = function (Tt)
        nt = nt_of(Tt)
        fc = 1 - (qt_sheath(Tt, nt) * Rt) / (qleg * Ru)
        return clamp(fc, ϵ, 1 - ϵ), nt
    end

    # R3 ⇒ g(T_t)=0 (minimize g^2)
    g(Tt) = begin
        fc, _ = fc_of(Tt)
        Tu^(3.5) - Tt^(3.5) - (7/2) * (Rt/Ru) * (qleg * fcond * Lpar) / (ke * (1 - fc))
    end

    function bracket_root(g, a, b; N::Int=200)
        xs = range(a, b; length=N)
        gprev = g(first(xs)); xprev = first(xs)
        for x in Iterators.drop(xs, 1)
            gx = g(x)
            if isfinite(gprev) && isfinite(gx) && sign(gprev) != sign(gx)
                return xprev, x
            end
            xprev = x; gprev = gx
        end
        return nothing
    end
    function bisect(g, a, b; tol=1e-8, maxit=100)
        fa = g(a); fb = g(b)
        if !(isfinite(fa) && isfinite(fb)) || fa*fb > 0
            return nothing
        end
        lo, hi, flo, fhi = a, b, fa, fb
        for _ in 1:maxit
            mid = 0.5*(lo+hi)
            fm  = g(mid)
            if !isfinite(fm); return nothing; end
            if abs(fm) < 1e-12 || (hi-lo) < tol
                return mid
            end
            if flo*fm <= 0
                hi, fhi = mid, fm
            else
                lo, flo = mid, fm
            end
        end
        return 0.5*(lo+hi)
    end
    Tt_lo = 0.1
    Tt_hi = 0.99*Tu

    ab = bracket_root(g, Tt_lo, Tt_hi; N=300)

    Tt_star = if ab !== nothing
        a, b = ab
        bisect(g, a, b; tol=1e-8, maxit=120)
    else
        xs = collect(range(Tt_lo, Tt_hi; length=400))
        vals = map(x->abs(g(x)), xs)
        xs[argmin(vals)]
    end

    fc_star, nt_star = fc_of(Tt_star)

    # write back & postprocess
    p.target.Te_t = Tt_star
    p.target.ne_t = nt_star
    p.coupling.fcool_optim = fc_star
    postprocess!(m)

    println("Direct solve (upstream_known) → Tt=$(Float64(Tt_star)), nt=$(Float64(nt_star)), fc=$(Float64(fc_star))")  

    return (minimizer = [Tt_star, nt_star/1e20, fc_star],
            minimum   = abs(g(Tt_star)),
            converged = (ab !== nothing))  
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