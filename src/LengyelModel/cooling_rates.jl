def Lint_ADAS(Tmin, Tmax,Texp,Lexp, Zimp; N = 100):
# Perform weighted cooling rate integral over specified temperature interval
# Inputs:  Tmin   minimum temperature for integral (eV)
#          Tmax   maximum temperature for integral (eV)
#          Texp   Exponent for temperature weighting
#          Lexp   Exponent for cooling rate weighting
#          Zimp   Z of impurity (Ar: 18; Kr: 36; Xe: 54)

   
    adas_data = load_adas_cooling_rates(Zimp)
        
    dT = LinRange(Tmin,Tmax,N)
    Lz = np.interp(Te,adasdata[Zimp]['temp'],adasdata[Zimp]['Lztot'])
    return np.trapz((Te^Texp)*(Lz^Lexp),Te) *1.e-6
end



get_cooling_rates(;database=:ADPAK)

get_adas_cooling_rates(imp::Symbol)




"""
    Calculates the appropriate neTau value to use for non-coronal adjustement to cooling rate.
    neTau increases with Te in a step-wise manner, with smoothing between steps.

    INPUTS:
    Te,  # array, Electron temperature [eV]

    OUTPUTS:
    neTau, # array, non-equilibrium parameter ne * tau [m^-3.s^-1]
    """
    using Interpolations
function Kallenbach_neTau(Te)
    Temodel = 10.0 .^ (collect(1:60)./10.0 .- 1)
    neTau_model = Temodel*0.
    neTau_model[Temodel.<=50] .= 1.e17
    neTau_model[(Temodel .> 50. .&& Temodel .< 100)] .= 1.e17 .+ (Temodel[(Temodel .> 50. .&& Temodel .< 100)] .- 50.0).*(1.e18-1.e17)/50.0
    neTau_model[Temodel .>=100 .&& Temodel .<=1000] .= 1.e18
    neTau_model[Temodel .>1000 .&& Temodel .< 2000] = 1.e18 .+ (Temodel[Temodel .>1000 .&& Temodel .< 2000] .- 1000.) .* (1.e19-1.e18)/1000.0
    neTau_model[Temodel .>= 2000 .&& Temodel .<= 100000] .= 1.e19
    interp_linear = Interpolations.linear_interpolation(Temodel, neTau_model)
    return interp_linear(Te)
end
""" 
Perform weighted cooling rate integral over specified temperature interval
# Inputs:  Tmin   minimum temperature for integral (eV)
#          Tmax   maximum temperature for integral (eV)
#          Texp   Exponent for temperature weighting
#          Lexp   Exponent for cooling rate weighting
#          Zimp   Z of impurity (Ar: 18; Kr: 36; Xe: 54)
"""

function Lint(Tmin, Tmax,Texp,Lexp, imp::Symbol; N=101, cooling_rate_data=get_cooling_rates())

    Te = Tedata
    Lz = Lzdata*1.e-6
    T = LinRange(Tmin,Tmax,N)
    Lz_ = Interpolations.linear_interpolation(Te,Lz)
    return NumericalIntegration.integrate(T,T .^ Texp .* Lz_(T).^ Lexp)
end
    


using Plots
fraction =collect(0.01:0.01:0.05)
plot()
for f in fraction
v = []
Tup = collect(0.0:10:200.0)
for Tu in Tup
    #V_legyel__ADAS(0.0, Tu, f, :Ne)
    push!(v,(V_legyel__ADAS(0.0, Tu, f, :Ne)/1.6e-19))
end

plot!(Tup, v; label="f=$f")
end
layout = @layout[
  p1 [
    grid(1, 1)
    p
  ]
]

plot(;layout= layout)

v = []
nu = 1e20
qtarget = 1e7
Tup = 200.0
fraction =collect(0.001:0.001:0.2)
for f in fraction
    push!(v,sqrt((nu*Tup*V_legyel__ADAS(0.0, Tup, f, :Ar))^2+qtarget^2)*1e-6*1.35e-3*(2*pi))
end

plot!(v,fraction; label="Ar" )
v = []
nu = 1e20
qtarget = 1e7
Tup = 200.0
fraction =collect(0.001:0.001:0.2)
for f in fraction
    push!(v,sqrt((nu*Tup*V_legyel__ADAS(0.0, Tup, f, :Ne))^2+qtarget^2)*1e-6*1.35e-3*(2*pi))
end
plot!(v,fraction; label="Ne" )
using Plots
pythonplot()
using Colors
colors = palette(:tab10)


plot()
imps =[:Ne,:Ar,:Kr,:Xe]
for (i,imp) in enumerate(imps)
v = []
v2 = []
zeffup = []
nu = 1e20
qtarget = 1e7
Tup = 200.0
fraction =10 .^ LinRange(-4:0.05:-1)

for f in fraction
    zeff = ADAS.get_effective_charge(imp);
    push!(zeffup,zeff(f,nu,Tup))
    push!(v,sqrt((nu*Tup*V_legyel__ADAS(0.0, Tup, f, imp))^2+qtarget^2)*1e-6*1.35e-3*(2*pi))
    push!(v2,sqrt((3*nu*Tup*V_legyel__ADAS(0.0, Tup, f, imp))^2+qtarget^2)*1e-6*1.35e-3*(2*pi))
end

p = plot!(v,zeffup; label="$imp ", linewidth = 3.0, color = colors[i])
plot!(v2,zeffup; label=nothing,linewidth = 3.0, linestyle=:dash, color = colors[i])
end
plot!([1], [0], label = "Lengyel model: GASC", color = "black")
plot!([1], [0], linestyle = :dash, label = "Lengyel model: 1D/2D informed", color = "black")
plot!(xlabel=L"$P_{SOL} B_θ / R_{omp}#$", 
ylabel = L"$Z_{eff}#$ at separatrix upstream", 
legend_position=:outertop, 
frame=true,xlim=[1,25], 
ylim=[1.0,5.0], 
aspect_ratio = (24)/(4),
xlabelfontsize=20, 
ylabelfontsize=20, 
size=(600,600), 
legendfontsize=12, 
legend_columns=1,
xtickfontsize=15,
ytickfontsize=15
)



imps =[:Ne,:Ar,:Kr,:Xe]
for imp in imps
v = []
zeffup = []
nu = 1e20
qtarget = 1e7
Tup = 200.0
fraction =collect(0.001:0.001:0.2)

for f in fraction
    zeff = ADAS.get_effective_charge(imp);
    push!(zeffup,zeff(f,nu,Tup))
    
end

plot!(v,zeffup; linestyle=:dash ,label="$imp [Moulton2021]" )
end