using SimulationParameters: AbstractParameters, Entry, Switch

Base.@kwdef mutable struct LengyelModelPlasmaParameters{T<:Real} <: AbstractParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    P_SOL::Entry{T} = Entry{T}("W", "Power coming through the separatrix")
    R_omp::Entry{T} = Entry{T}("m", "Outer midplane radius")
    Ip::Entry{T} = Entry{T}("A", "Plasma current")
    κ::Entry{T} = Entry{T}("-", "Elongation")
    ϵ::Entry{T} = Entry{T}("-", "Plasma aspect ratio")
    Bt_omp::Entry{T} = Entry{T}("T", "Toroidal magnetic field at the outer midplane")
    Bpol_omp::Entry{T} = Entry{T}("T", "Poloidal magnetic field at the outer midplane")
end

Base.@kwdef mutable struct LengyelModelSOLParameters{T<:Real} <: AbstractParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    n_up::Entry{T} = Entry{T}("m⁻³", "Upstream density")
    T_up::Entry{T} = Entry{T}("eV", "Upstream temperature")
    λ_omp::Entry{T} = Entry{T}("m", "Heat flux decay length")
    f_imp::Entry{Vector{Float64}} = Entry{Vector{Float64}}("-", "Vector of impurity fractions")
    imp::Entry{Vector{Symbol}} = Entry{Vector{Symbol}}("-", "Vector of impurity species")
    f_adhoc::Entry{T} = Entry{T}("-", "Radiation enhancement factor"; default=1.0)
end

Base.@kwdef mutable struct LengyelModelTargetParameters{T<:Real} <: AbstractParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    f_omp2target_expansion::Entry{T} = Entry{T}("-", "Flux expansion and projection of λ_omp onto the target")
    R::Entry{T} = Entry{T}("m", "Major radius at the outer strike point")
    f_spread_pfr::Entry{T} = Entry{T}("-", "Heat flux expansion factor in the private flux region (eg. due to transport) should be >= 1.0")
    α_sp::Entry{T} = Entry{T}("rad", "Pitch angle at the outer strike point")
    θ_sp::Entry{T} = Entry{T}("rad", "Poloidal angle of the target at the outer strike point") # [CURRENTLY NOT USED]
end

Base.@kwdef mutable struct LengyelIntegralParameters{T<:Real} <: AbstractParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    T_down::Entry{T} = Entry{T}("eV", "Temperature downstream"; default=0.0)
    Zeff_exp::Entry{T} = Entry{T}("-", "Exponent of the Zeff dependency in Spitzer heat conductivity"; default=-0.3)
    Texp::Entry{T} = Entry{T}("-", "Exponent of the temperature dependency in Spitzer heat conductivity"; default=0.5)
    Lexp::Entry{T} = Entry{T}("-", "Exponent of the cooling dependency in Spitzer heat conductivity"; default=1.0)
    κ0::Entry{T} = Entry{T}("W⋅m/eV^3.5", "Heat conductivity"; default=2390.0)
end


Base.@kwdef mutable struct LengyelModelParameters{T<:Real} <:DivertorHeatFluxModelParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :LengyelModelParameters
    sol::LengyelModelSOLParameters{T} = LengyelModelSOLParameters{T}()
    plasma::LengyelModelPlasmaParameters{T} = LengyelModelPlasmaParameters{T}()
    target::LengyelModelTargetParameters{T} = LengyelModelTargetParameters{T}()
    integral::LengyelIntegralParameters{T} = LengyelIntegralParameters{T}()
end

Base.convert(::Type{DivertorHeatFluxModelParameters{Float64}}, t::T) where {T<:DivertorHeatFluxModelParameters} = t