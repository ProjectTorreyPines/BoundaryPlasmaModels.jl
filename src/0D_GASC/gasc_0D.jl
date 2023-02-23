import numpy as np
from gasc.cooling_models import ( get_cooling_rates, load_adas_cooling_rates, read_adas_file, Kallenbach_neTau, \
									 Mavrin_cooling_rates, Lint, Lint_Mavrin, Lint_ADAS )
μ₀ = 4.e-7 * pi

compute_q⋆(R::T,κ::T,Ip::T,ϵ::T,Bt::T) where {T<:Real} = pi*(ϵ*R)^2. * (1. + κ^2.) * Bt/(μ₀*R*Ip*1.e6)
compute_ℒ_omp_target(R::T,κ::T,Ip::T,ϵ::T,Bt::T) where {T<:Real} = 4.33 * compute_q⋆(R,κ,Ip,ϵ,Bt) * R
compute_A⟂(R::T,ϵ::T,λ_q::T, Bpol::T, Bt::T) where {T<:Real} = 2.*pi*R*(1.+epsilon)*λ_q *sin(arctan(Bpol/Bt))
compute_q_∥∥_omp(P_SOL::T, R::T,ϵ::T,λ_q::T, Bpol::T, Bt::T) = P_SOL*1.e6/compute_A⟂(R,ϵ,λ_q, Bpol, Bt)

function λq_sieglin(R,P_SOL,Bpol,ϵ,νₑ)
    """Eich scaling(NF 53 093031) & B. Sieglin PPCF 55 (2013) 124039"""
    # Eich scaling(NF 53 093031)
    λ_q = 1.35 * (P_SOL)^(-0.02) * R^0.04 * Bpol^(-0.92) * νₑ^0.42 * 1.e-3

    #From B. Sieglin PPCF 55 (2013) 124039
    S = 0.09*(νₑ/1.e19)^1.02 * Bpol^(-1.01)*(R/3.0) * 1.e-3

    return λ_q + 1.64*S
end
   """Eich scaling (NF 53 093031)"""
λq_eich(R::T,P_SOL::T,Bpol::T,ϵ::T) = 1.35e-3 * P_SOL^(-0.02)* R^0.04 * Bpol^(-0.92) * ϵ^(0.42) # Eich scaling (NF 53 093031)

λq_loarte(P_SOL::T,B0::T,q95iter::T)  = 0.00265 * P_SOL^0.38 * B0^(-0.71) * q95iter ^ 0.3

def model_equations (x, *args):
    # xorder:  Tup,qI,Lr,Tsh,qsh,fm,qpl
    # args_order:  qup, fz, nup, L_parallel
    
    # collect passed info into variables
    Tmid = x[0]**2.
    qI = x[1]**2.
    Lr = x[2]**2.
    Tsh = x[3]**2.
    qsh = x[4]**2.
    #fm = x[5]**2.
    #qpl = x[6]**2.
    #Tsh = x[1]**2.
    #Lm = x[2]**2.
    #qplate = x[3]**2.

    qmid = args[0]
    fz = args[1]
    nmid = args[2]
    Lparallel = args[3]
    Zimp = args[4]

    f = x*0.
    frad_correction = 1.0

    ee = 1.602e-19
    mi = 1.67e-27
    cs0 = np.sqrt(2.*ee/mi)
    Tc = np.maximum(Tsh,15.)
    #Tc = 450.
    #Tsh = np.maximum(Tsh,Tc)
    gamma = 9.2
    etaB = 1.0
    kappa0 = 2390.
    xi = 2.8 * np.exp(-13.6/Tsh)*Tsh**0.19/(6.+0.0873*Tsh)
    alpha = xi / (1. + xi)
    #fm = 1. - 2.*(alpha/(1.+alpha))**((alpha+1)/2.)
    #qpl = (1.-fm)*gamma/2./etaB*ee*nmid*Tmid*cs0*np.sqrt(Tsh)
    Lm = Lparallel-Lr

    f[0] = qmid**2. - qI**2. - 2.*kappa0*(nmid*Tmid)**2.*fz*Lint(Tc,Tmid,0.5,1.0,Zimp)*frad_correction
    f[0] = f[0]

    dT = (Tmid-Tc)/101
    nT = 101
    Teint = Tc + np.arange(nT)*dT
    qT = np.zeros(nT)
    for i in range(nT):
       qT[i] = np.sqrt(qI**2. + 2.*kappa0*(nmid*Tmid)**2.*fz*Lint(Tc,Teint[i],0.5,1.0,Zimp)*frad_correction)

    intgrnd = (Teint**2.5)/qT
    f[1] = Lr - kappa0 * np.trapz(intgrnd,x=Teint)

    #f[2] = qI - gamma/2.*ee*nmid*Tmid*cs0*np.sqrt(Tc)
    #nT = 101
    #dT = (Tc-Tsh)/nT
    #Teint = Tsh + np.arange(nT)*dT
    #intgrnd = Teint**1.5/np.interp(Teint,Te,Lz)
    #if (Tc-Tsh) < 0.0:
    #    intgrnd = Teint*0.


    f[2] = qI - gamma/2.*ee*nmid*Tmid*cs0*np.sqrt(Tc)
    f[3] = Lint(Tsh,Tc,1.5,-1.0,Zimp) - (nmid*Tmid*Lm*fz)/ (gamma*ee*cs0)
    f[4] = qsh - gamma/2.*ee*nmid*Tmid*cs0*np.sqrt(Tsh)
    #f[5] = fm - (1. - 2*(alpha/(1.+alpha))**((alpha+1)/2.))
    fm = 1. - 2*(alpha/(1.+alpha))**((alpha+1)/2.)

    #f[6] = qpl - (1.-fm)*qsh/etaB
    qpl = (1.-fm)*qsh/etaB

    return f

def SOL_model_equations (x, *args):
#xorder:  Tmid, qplate
#args_order:  qup, nup, L_parallel, Zeff, fz1,Zimp1, fz2,Zimp2

    Tmid = x[0]
    qplate = x[1]

    qmid = args[0]
    nmid = args[1]
    Lparallel = args[2]
    effectiveZ = args[3]
    fz1 = args[4]
    Zimp1 = args[5]
    fz2 = args[6]
    Zimp2 = args[7]
    f = x*0.

    Tc = 0.
    kappa0 = 2390.*effectiveZ**(-0.3)
    f[0] = qmid**2. - qplate**2. - 2.*kappa0*(nmid*Tmid)**2.*(fz1*Lint(Tc,Tmid,0.5,1.0,Zimp1)+fz2*Lint(Tc,Tmid,0.5,1.0,Zimp2))

    dT = (Tmid-Tc)/101
    nT = 101
    Teint = Tc + np.arange(nT)*dT
    qT = np.zeros(nT)
    for i in range(nT):
       qT[i] = np.sqrt(qmid**2. - 2.*kappa0*(nmid*Tmid)**2.*(fz1*Lint(Tc,Teint[i],0.5,1.0,Zimp1)+fz2*Lint(Tc,Tmid,0.5,1.0,Zimp2)))

    intgrnd = kappa0*Teint**2.5/qT
    f[1] = Lparallel - np.trapz(intgrnd,Teint)

    return f

def broyfunc (x, *args):

    qu = args[0]
    nu = args[1]
    fzSol1 = args[2]
    Zimp1 = args[3]
    fzSol2 = args[4]
    Zimp2 = args[5]
    Lparallel = args[6]
    m_D = args[7]
    eps_pot = args[8]
    eps_rad = args[9]
    gamma = args[10]
    kappa _0 = args[11]

    #m_D =
    #Lz = 2.5e-33
    #kappa_0 = 1.3e69 * 1.602e-19**3.5  * 2**(-0.3)
    #eps_pot = 15.
    #eps_rad = 16.
    #gamma = 7.
    #L = args[3]
    Tu = x[0]**2.
    Tl = x[1]**2.
    nt = x[2]**2.

    Tgrid = np.arange(101)*(Tu-Tl)/100. + Tl
    neTau = Kallenbach_neTau(Tgrid)
    neTau = neTau*0. + 1.e19
    #Lz = np.trapz(Mavrin_cooling_rates(Tgrid,neTau,18)*np.sqrt(Tgrid),Tgrid)/np.trapz(np.sqrt(Tgrid),Tgrid)

    #try:
    #    frad = 1.-np.sqrt(1.-14./3.*nu**2.*Cz*Lz*L/qu)
    #except:
    #    frad = 1.

    #qtold = qu * (1.-frad)
    tmp = qu**2. - 2.*kappa_0*nu**2.*Tu**2.*(fzSol1*np.trapz(Mavrin_cooling_rates(Tgrid,neTau,Zimp1)*np.sqrt(Tgrid),Tgrid)+
                                             fzSol2*np.trapz(Mavrin_cooling_rates(Tgrid,neTau,Zimp2)*np.sqrt(Tgrid),Tgrid))
    
    q_rad = np.sqrt(2.*kappa_0*nu**2.*Tu**2.*(fzSol1*np.trapz(Mavrin_cooling_rates(Tgrid,neTau,Zimp1)*np.sqrt(Tgrid),Tgrid)+
                                              fzSol2*np.trapz(Mavrin_cooling_rates(Tgrid,neTau,Zimp2)*np.sqrt(Tgrid),Tgrid)))

    if tmp >  0:
        qt = np.sqrt(tmp)
    else:
        qt = 0.

    f1 = Tu**3.5 - Tl**3.5 - 7.*qu*Lparallel/(2.*kappa_0)
    Tsh = Tl

    xi = 2.8 * np.exp(-13.6/Tsh)*Tsh**0.19/(6.+0.0873*Tsh)
    alpha = xi / (1. + xi)
    fm = 2*(alpha/(1.+alpha))**((alpha+1)/2.)
    f2 = fm*nu*Tu - 2.*nt*Tl
    f3 = qt - nt * np.sqrt(2.*Tl*1.602e-19/m_D)*(gamma*Tl + eps_pot + eps_rad)*1.602e-19

    return ([f1,f2,f3])

        
struct GASC0DHeatFluxParameters
    P_SOL :: Float64 # power into SOL [MW]
    R :: Float64 # major radius [m]
    Ip         # plasma current [MA]
    κ      # elongation
    ϵ    # inverse aspect ratio
    Bt   # toroidal field to use for flux area [T]
    Bpol  # poloidal field to use for flux area [T]
    Ndiv       # number of divertors
end

kappa,      # elongation

end
def qparallel_divertor(
    Psol,       
    lambda_int, # heat flux width to use for flux area [m]
    nu,         # upstream density []
    fzHe,       # Helium concentration in divertor
    ZHe,        # Helium atomic number
    fz1Sol,     # impurity 1 concentration in divertor
    Zimp1,      # impurity 1 atomic number
    fz2Sol,     # impurity 2 concentration in divertor
    Zimp2,      # impurity 2 atomic number
    Zeff, # effective Z in divertor
    ):

    '''
    Follows the two-point model (Kallenbach PPCF 55 (2013) 124041) to calculate the
    heat flux at the divertor target with impurity radiation.

    Psol,       # power into SOL [MW]
    R,          # major radius [m]
    Bt,         # toroidal field to use for flux area [T]
    Bpol,       # poloidal field to use for flux area [T]
    lambda_int, # heat flux width to use for flux area [m]
    Ip,         # plasma current [MA]
    epsilon,    # inverse aspect ratio
    kappa,      # elongation
    Ndiv,       # number of divertors
    nu,         # upstream density []
    fzHe,       # Helium concentration in divertor
    ZHe,        # Helium atomic number
    fz1Sol,     # impurity 1 concentration in divertor
    Zimp1,      # impurity 1 atomic number
    fz2Sol,     # impurity 2 concentration in divertor
    Zimp2,      # impurity 2 atomic number
    effectiveZ, # effective Z in divertor
    '''

    #_lambdaq = 1.35 * (Psol)**(-0.02)*R**0.04*Bpol**(-0.92)*epsilon**0.42 *1.e-3

    #From B. Sieglin PPCF 55 (2013) 124039
    #S = (0.09*(nu/1.e19)**1.02*Bpol**-1.01)*(R/3.0) * 1.e-3

    #lambda_int = _lambdaq + 1.64*S

    # Calculate q_parallel at the outboard midplane (Kallenback Eqn. 1)

    # Calculate connection length from outboard midplane to divertor target

    # Calculate upstream temperature 
    ℒ_∥∥ = compute_ℒ_omp_target(R,κ,Ip,ϵ,Bt)
    q_∥∥ = compute_q_∥∥_omp(R,ϵ,λ_q, Bpol, Bt)
    κ₀ = 2390.* Z_eff^(-0.3) 
    compute_Tu( q_∥∥, ℒ_∥∥,κ₀) = (7.*q_∥∥*ℒ_∥∥/(2.*κ₀))^(2./7.) 

    # Integrate to get power radiated along flux tube (Kallenbach Eqn. 5) 
    # Lint_ADAS(Tmin,Tmax,Texp,Lexp,Zimp) performs the weighted cooling rate integral over the specified
    # temperature interval. 
    q_rad =  nu*Tu*np.sqrt(2.*kappa_0*(fzHe*Lint_ADAS(0.,Tu,0.5,1.0,ZHe)+
                                       fz1Sol*Lint_ADAS(0.,Tu,0.5,1.0,Zimp1)+
                                       fz2Sol*Lint_ADAS(0.,Tu,0.5,1.0,Zimp2) ))

    #Tgrid = np.arange(101)*Tu/101.

    #Lz1 = Lint(0.,Tu,0.5,1.0,Zimp1)/Lint(0.,Tu,0.5,0.0,Zimp1)
    #Lz2 = Lint(0.,Tu,0.5,1.0,Zimp2)/Lint(0.,Tu,0.5,0.0,Zimp2)

    #frad = 1.-np.sqrt(1.- 14./3.*nu**2.*(fz1Sol*Lz1+fz2Sol*Lz2)*L_parallel/q_parallel)

    #q_rad2 =  np.sqrt(2.*kappa_0*effectiveZ**(-0.3)*(nu*Tu)**2.*(fz1Sol*Lint(0.,Tu,0.5,1.0,Zimp1)+fz2Sol*Lint(0.,Tu,0.5,1.0,Zimp2) ))

    q_target = np.sqrt(q_parallel**2. - q_rad**2.)   #Expected heat flux without geometric considerations

    P_rad = q_rad * A_perp * Ndiv * 1.e-6

    return (q_parallel, q_rad, q_target, Tu, P_rad) #, lambda_int)


    m_D = 1.67e-27
    kappa_0 = 2390.  * effectiveZ**(-0.3)
    eps_pot = 15.
    eps_rad = 16.
    gamma = 7.

    args = [q_parallel,nu,fz1Sol,Zimp1,fz2Sol,Zimp2,L_parallel,m_D,eps_pot,eps_rad,gamma,kappa_0]
    xguess = np.array([Tu,10.,nu])
    solution = sc.root(broyfunc,np.sqrt(xguess),args=tuple(args),tol=1.e-6,method='lm')
    if solution['success']:
        Tu = solution['x'][0]**2.
        Tt = solution['x'][1]**2
        nt = solution['x'][2]**2.
        q_target_old = q_target
        q_target = nt * np.sqrt(2.*Tt*1.602e-19/m_D)*(gamma*Tt + eps_pot + eps_rad)*1.602e-19
        #q_rad_old = q_rad
        q_rad = np.sqrt(q_parallel**2 - q_target**2.)
        P_rad = q_rad * A_perp * Ndiv * 1.e-6

    f = model_equations(np.sqrt(xguess),*args)

    solution = sc.fsolve(model_equations,np.sqrt(xguess),args=tuple(args),xtol=1.e-6)
    q_target = solution[1]
    q_rad = np.sqrt(q_parallel**2.-q_target**2.)
    Tu = solution[0]

    return (q_parallel, q_rad, q_target, Tu, P_rad, lambda_int)


function compute fz_Reinke(P_SOL, R, Bt, Bpol, Ip, epsilon, κ , nu, Zimp):

    λ_q = 1.35 * P_SOL^(-0.02) *R ^0.04* Bpol^(-0.92) * ϵ^0.42 * 1.e-3
    q_∥∥ = P_SOL*1.e6 * sqrt(Bt^2.+Bpol^2.) / (2.*pi*R*λ_q*Bpol)
    #q_parallel2 = Psol*1.e6/(2.*np.pi*R*lambdaq*np.sin(np.arctan(Bpol/Bt)))

    q⋆ = compute_q⋆(R,κ,Ip,ϵ,Bt)
    ℒ_∥∥ = 4.33*q⋆*R

    #q_parallel2 = 117.9*Psol*np.sqrt(Bt**2.+Bpol**2.)/ (R * epsilon**0.42)*1.e6

    κ₀ = 1.3e69 * (1.602e-19)^3.5
    #kappa_0 = 3.09e3 / Zimp /14.  #mks units with T in eV (Stacey Fusion Plasma Physics pg. 379)  - Using 14. for ln(gamma)
      #eV
    Tu = compute_Tu( q_∥∥, ℒ_∥∥,κ₀)    #eV

    fz = q_∥∥^2./(2.*κ₀*nu^2.*Tu^2.*Lint(0.,Tu,0.5,1.0,Zimp))

    return fz