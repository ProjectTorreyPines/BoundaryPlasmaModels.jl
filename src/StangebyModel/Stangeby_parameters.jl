using SimulationParameters: AbstractParameters, Entry, Switch

Base.@kwdef mutable struct StangebyModelPlasmaParameters{T<:Real} <: AbstractParameters{T}
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

Base.@kwdef mutable struct StangebyModelSOLParameters{T<:Real} <: AbstractParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    T_up::Entry{T} = Entry{T}("eV", "Upstream temperature")
    λ_omp::Entry{T} = Entry{T}("m", "Heat flux decay length")
end

Base.@kwdef mutable struct StangebyModelTargetParameters{T<:Real} <: AbstractParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    L_para::Entry{T} = Entry{T}("m", "Connection length")
    f_omp2target_expansion::Entry{T} = Entry{T}("-", "Flux expansion and projection of λ_omp onto the target")
    R::Entry{T} = Entry{T}("m", "Major radius at the outer strike point")
    f_spread_pfr::Entry{T} = Entry{T}("-", "Heat flux expansion factor in the private flux region (eg. due to transport) should be >= 1.0")
    α_sp::Entry{T} = Entry{T}("rad", "Pitch angle at the outer strike point")
    θ_sp::Entry{T} = Entry{T}("rad", "Poloidal angle of the target at the outer strike point") # [CURRENTLY NOT USED]
end

Base.@kwdef mutable struct StangebyIntegralParameters{T<:Real} <: AbstractParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    f_cond::Entry{T} = Entry{T}("-", "fraction of conducted heat flux"; default=1.0)
    κ0::Entry{T} = Entry{T}("W⋅m/eV^3.5", "Heat conductivity"; default=2390.0)
end

Base.@kwdef mutable struct StangebyModelParameters{T<:Real} <: DivertorHeatFluxModelParameters{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :StangebyModelParameters
    sol::StangebyModelSOLParameters{T} = StangebyModelSOLParameters{T}()
    plasma::StangebyModelPlasmaParameters{T} = StangebyModelPlasmaParameters{T}()
    target::StangebyModelTargetParameters{T} = StangebyModelTargetParameters{T}()
    integral::StangebyIntegralParameters{T} = StangebyIntegralParameters{T}()
end
