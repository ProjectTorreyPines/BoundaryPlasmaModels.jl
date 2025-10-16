import ..BoundaryPlasmaModels: SOLBoundaryParameters
using SimulationParameters: AbstractParameters, Entry, Switch

Base.@kwdef mutable struct SLCoupledPlasmaParameters{T<:Real} <: AbstractParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    P_SOL::Entry{T}   = Entry{T}("W", "Power coming through the separatrix")
    f_psol::Entry{T}  = Entry{T}("-",  "Fraction of Psol"; default=1.0)
    R_omp::Entry{T}   = Entry{T}("m", "Outer midplane radius")
    R_x::Entry{T}     = Entry{T}("m", "X-point radius")
    Ip::Entry{T}      = Entry{T}("A", "Plasma current")
    κ::Entry{T}       = Entry{T}("-", "Elongation")
    ϵ::Entry{T}       = Entry{T}("-", "Plasma aspect ratio")
    Bt_omp::Entry{T}  = Entry{T}("T", "Toroidal magnetic field at the outer midplane")
    Bpol_omp::Entry{T}= Entry{T}("T", "Poloidal magnetic field at the outer midplane")
end

Base.@kwdef mutable struct SLCoupledSOLParameters{T<:Real} <: AbstractParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :SLCoupledSOLParameters
    λ_omp::Entry{T}  = Entry{T}("m",   "Heat flux decay length at OMP")
    Te_up::Entry{T}  = Entry{T}("eV",  "Upstream electron temperature (solved)")
    ne_up::Entry{T}  = Entry{T}("m^-3","Upstream electron density (solved)")
end

Base.@kwdef mutable struct SLCoupledTargetParameters{T<:Real} <: AbstractParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :SLCoupledTargetParameters
    L_para::Entry{T}                 = Entry{T}("m",   "Parallel connection length")
    f_omp2target_expansion::Entry{T} = Entry{T}("-",   "Flux expansion OMP→target")
    R_strike::Entry{T}               = Entry{T}("m",   "Geometric major radius at strike point")
    f_spread_pfr::Entry{T}           = Entry{T}("-",   "Heat-flux spread factor in PFR (>=1)"; default = 1.0)
    α_sp::Entry{T}                   = Entry{T}("rad", "Toroidal tilt angle at strike point")
    θ_sp::Entry{T}                   = Entry{T}("rad", "Poloidal tilt angle at strike point")
    Te_t::Entry{T}                   = Entry{T}("eV",  "Target electron temperature (if available)"; default = 5.0)
    ne_t::Entry{T}                   = Entry{T}("m^-3","Target electron density (if available)";     default = 1e19)
end

Base.@kwdef mutable struct SLCoupledCouplingParameters{T<:Real} <: AbstractParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :SLCoupledCouplingParameters
    f_cond::Entry{T}     = Entry{T}("-",  "Fraction of conductive heat"; default=1.0)
    kappa_e::Entry{T}    = Entry{T}("W⋅m/eV^3.5", "Effective e- conduction prefactor"; default=2390.0)
    mi_amu::Entry{T}     = Entry{T}("amu","Main ion mass"; default=2.0)
    tau_t::Entry{T}      = Entry{T}("-",  "Ti/Te at target"; default=1.0)
    tau_u::Entry{T}      = Entry{T}("-",  "Ti/Te upstream"; default=1.0)
    M_t::Entry{T}        = Entry{T}("-",  "Mach at target"; default=1.0)
    M_u::Entry{T}        = Entry{T}("-",  "Mach upstream"; default=0.0)
    # log10(1 - fcool) = a0 + a1 log10 Te_t + ... + a4 (log10 Te_t)^4
    fcool_optim::Entry{T}          = Entry{T}("-",  "optimized fcool"; default=0.0)
    fcool_poly::Entry{NTuple{5,T}} = Entry{NTuple{5,T}}("-", "Poly for log10(1 - f_cool) vs log10 Te_t"; default=(-1.171, 2.101, -3.212, 3.567, -1.477))
    fmom_poly::Entry{NTuple{5,T}}  = Entry{NTuple{5,T}}("-", "Poly for log10(1 - f_mom) vs log10 Te_t"; default=(-0.333, 1.173, -2.020, 1.972, -0.751))
end

Base.@kwdef mutable struct SLCoupledRadiationParameters{T<:Real} <: AbstractParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :SLCoupledRadiationParameters
    species::Entry{Vector{Symbol}} = Entry{Vector{Symbol}}("-", "Impurity species"; default=Symbol[])
    fractions::Entry{Vector{T}}    = Entry{Vector{T}}("-", "Impurity fractions"; default=T[])
    impurity_fraction::Entry{Vector{T}}   = Entry{Vector{T}}("-", "impurity fraction ( n_z / n_e)"; default=T[])
    q_rad::Entry{T}   = Entry{T}("W/m²", "Radiated heat flux density (solved)"; default= NaN)
    cz_values::Entry{Vector{T}}    = Entry{Vector{T}}("-", "cz per species (solved)"; default=T[])
    #Zeff::Entry{T}    = Entry{T}("-",   "Effective charge (solved)"; default=NaN)
end

Base.@kwdef mutable struct SLCoupledParameters{T<:Real} <: SOLBoundaryParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :SLCoupledParameters
    plasma::SLCoupledPlasmaParameters{T}      = SLCoupledPlasmaParameters{T}()
    sol::SLCoupledSOLParameters{T}            = SLCoupledSOLParameters{T}()
    target::SLCoupledTargetParameters{T}      = SLCoupledTargetParameters{T}()
    coupling::SLCoupledCouplingParameters{T}  = SLCoupledCouplingParameters{T}()
    radiation::SLCoupledRadiationParameters{T}= SLCoupledRadiationParameters{T}()
end
