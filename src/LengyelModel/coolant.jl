
import Base:Iterators
using Plots
using LaTeXStrings
gr()
#import Plots.PythonPlot.pyplot as pyplt
function meshgrid(x,y)
    return first.(Iterators.product(x,y)), last.(Iterators.product(x,y))
end
clims = (800.0,2000.0)
#PythonPlot.figure()

# ITER Baseline water 
# paper: 
# [Bonin] Bonnin, X., Pitts, R. A., Komarov, V., Escourbiac, F., Merola, M., Bo, L., ... & Kukushkin, A. S. (2017). ITER divertor plasma response to time-dependent impurity injection. Nuclear Materials and Energy, 12, 1100-1105. 
# [Morshedy] El-Morshedy, Salah El-Din. "Thermal-hydraulic modelling and analysis of ITER tungsten divertor monoblock." Nuclear Materials and Energy 28 (2021): 101035.

δ = collect(0.001:0.001:0.01)
T_inlet = 273 + 150
α = collect(1:0.01:9.0).*1e4
α_,δ_ =  meshgrid(α,δ)
κ = 164
q = 1e7
f_circular = 1.07
T_max = @. (q/κ*δ_ + q/α_ + T_inlet) * f_circular

p1= contourf(δ,α, T_max;colorbar=false, clim = clims,)
cs = contour!(δ,α, T_max; levels=1200:1200, linewidth=3, contour_labels = true, label="1200K",fontsize=10, color = cgrad([:red,:red,:red]))
cs = contour!(δ,α, T_max; levels=1000:1000, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="1000K")
cs = contour!(δ,α, T_max; levels=1500:1500, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="1500K")
cs = contour!(δ,α, T_max; levels=2000:2000, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="2000K")
plot!(;xlabel=L"\delta_{surface \rightarrow pipe} [m]")
plot!(;ylabel=L" α $[W.m^{-2}.K^{-1}]$")
vline!([0.006];color=:green, linewidth=2,label="ITER monoblock design [Bonin]")
vline!([0.008];color=:green, linewidth=2,label="ITER monoblock design [Morshedy]")

hline!([0.08e6];color=:green, linewidth=2,linestyle = :dash, label="water at Q= 10 MW.m^{-2}[Morshedy]")
plot!(;legend_position=legend_position=:outertop)
plot!(title=L"ITER cooling: Q_{target} = %$(q/1e6) MW; $T_{inlet}$ = %$(T_inlet) K, $f_{circular}$ = %$(f_circular)")
plot!(titlefontsize=10)


#%%
δ = collect(0.001:0.001:0.01)
T_inlet = 273 + 150
α = collect(1:0.01:9.0).*1e4
α_,δ_ =  meshgrid(α,δ)
κ = 164
q = 2e7
f_circular = 1.07
T_max = @. (q/κ*δ_ + q/α_ + T_inlet) * f_circular

p2= contourf(δ,α, T_max;colorbar=false, clim = clims)
cs = contour!(δ,α, T_max; levels=1200:1200, linewidth=3, contour_labels = true, label="1200K",fontsize=10, color = cgrad([:red,:red,:red]))
cs = contour!(δ,α, T_max; levels=1000:1000, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="1000K")
cs = contour!(δ,α, T_max; levels=1500:1500, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="1500K")
cs = contour!(δ,α, T_max; levels=2000:2000, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="2000K")
plot!(;xlabel=L"\delta_{surface \rightarrow pipe} [m]")
plot!(;ylabel=L" α $[W.m^{-2}.K^{-1}]$")
vline!([0.006];color=:green, linewidth=2,label="ITER monoblock design [Bonin]")
vline!([0.008];color=:green, linewidth=2,label="ITER monoblock design [Morshedy]")

hline!([0.08e6];color=:green, linewidth=2,linestyle = :dash, label="water at Q= 10 MW.m^{-2}[Morshedy]")
plot!(;legend_position=legend_position=:outertop)
plot!(title=L"ITER cooling: Q_{target} = %$(q/1e6) MW; $T_{inlet}$ = %$(T_inlet) K, $f_{circular}$ = %$(f_circular)")
plot!(titlefontsize=10)

##%

# Helium coolant
# paper: 
# [Baxi] 

δ = collect(0.001:0.001:0.01)
T_inlet = 273 + 200
α = collect(1:0.01:9.0).*1e4
α_,δ_ =  meshgrid(α,δ)
κ = 164
q = 1e7
f_circular = 1.07
T_max = @. (q/κ*δ_ + q/α_ + T_inlet) * f_circular

p3= contourf(δ,α, T_max;colorbar=false, clim=clims)
cs = contour!(δ,α, T_max; levels=1200:1200, linewidth=3, contour_labels = true, label="1200K",fontsize=10, color = cgrad([:red,:red,:red]))
cs = contour!(δ,α, T_max; levels=1000:1000, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="1000K")
cs = contour!(δ,α, T_max; levels=1500:1500, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="1500K")
cs = contour!(δ,α, T_max; levels=2000:2000, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="2000K")
plot!(;xlabel=L"\delta_{surface \rightarrow pipe} [m]")
plot!(;ylabel=L" α $[W.m^{-2}.K^{-1}]$")
#vline!([0.006];color=:green, linewidth=2,label="ITER monoblock design [Bonin]")
vline!([0.008];color=:green, linewidth=2,label="ITER monoblock design [Morshedy]")
hline!([0.026e6];label="He: smooth Cu tube[Baxi]")
hline!([0.04e6];label="He: jet Cu tube [Baxi]")
hline!([0.05e6];label="He: jet flow []")
plot!(;legend_position=legend_position=:outertop)
plot!(title=L"He cooling: Q_{target} = %$(q/1e6) MW; $T_{inlet}$ = %$(T_inlet) K, $f_{circular}$ = %$(f_circular)")

plot!(titlefontsize=10)


##%

δ = collect(0.001:0.001:0.01)
T_inlet = 273 + 200
α = collect(1:0.01:9.0).*1e4
α_,δ_ =  meshgrid(α,δ)
κ = 164
q = 2e7
f_circular = 1.07
T_max = @. (q/κ*δ_ + q/α_ + T_inlet) * f_circular
p4= contourf(δ,α, T_max;colorbar=false, clim=clims)
cs = contour!(δ,α, T_max; levels=1200:1200, linewidth=3, contour_labels = true, label="1200K",fontsize=10, color = cgrad([:red,:red,:red]))
cs = contour!(δ,α, T_max; levels=1000:1000, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="1000K")
cs = contour!(δ,α, T_max; levels=1500:1500, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="1500K")
cs = contour!(δ,α, T_max; levels=2000:2000, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="2000K")
plot!(;xlabel=L"\delta_{surface \rightarrow pipe} [m]")
plot!(;ylabel=L" α $[W.m^{-2}.K^{-1}]$")
#vline!([0.006];color=:green, linewidth=2,label="ITER monoblock design [Bonin]")
vline!([0.008];color=:green, linewidth=2,label="ITER monoblock design [Morshedy]")
hline!([0.026e6];label="He: smooth Cu tube[Baxi]")
hline!([0.04e6];label="He: jet Cu tube [Baxi]")
hline!([0.05e6];label="He: jet flow []")
plot!(;legend_position=legend_position=:outertop)
plot!(title=L"He cooling: Q_{target} = %$(q/1e6) MW; $T_{inlet}$ = %$(T_inlet) K, $f_{circular}$ = %$(f_circular)")

plot!(titlefontsize=10)
h2 = scatter([0,0], [0,1], zcolor=[0,4], clims=clims,
                 xlims=(1,1.1), label="", colorbar_title="Surface temperature [K]", framestyle=:none)

    l = @layout [grid(2, 2) a{0.1w}]
    p_all = plot(p1, p2, p3, p4, h2, layout=l, link=:all, size=(1300,1300))
    savefig(p_all, "~/test.png")
#p_all = plot(p1, p2, p3, p4, layout=(2, 2), link=:all)


plot()

δ = collect(0.001:0.001:0.01)
T_inlet = 273 + 200
α = collect(1:0.01:9.0).*1e4
α_,δ_ =  meshgrid(α,δ)
κ = 164
q = 1e7
f_circular = 1.07
T_max = @. (q/κ*δ_ + q/α_ + T_inlet) * f_circular

p3= contourf(δ,α, T_max;colorbar=false, clim=clims)
# cs = contour!(δ,α, T_max; levels=1200:1200, linewidth=3, contour_labels = true, label="1200K",fontsize=10, color = cgrad([:red,:red,:red]))
cs = contour!(δ,α, T_max; levels=1063:1063, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="ITER baseline: 1063K")
# cs = contour!(δ,α, T_max; levels=1500:1500, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="1500K")
# cs = contour!(δ,α, T_max; levels=2000:2000, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="2000K")
plot!(;xlabel=L"\delta_{surface \rightarrow pipe} [m]")
plot!(;ylabel=L" heat transfert coeff. $[W.m^{-2}.K^{-1}]$")
#vline!([0.006];color=:green, linewidth=2,label="ITER monoblock design [Bonin]")
vline!([0.006];color=:blue,linestyle=:dash, linewidth=2,label="ITER monoblock design [Morshedy]")
hline!([0.026e6];linewidth=2,label="He: smooth Cu tube[Baxi]")
# hline!([0.04e6];linewidth=2,label="He: jet Cu tube [Baxi]")
hline!([0.05e6];linewidth=2,label="He: jet flow []")
hline!([0.08e6];color=:blue, linestyle=:dot, linewidth=2,label="Water: ITER")
plot!(;legend_position=legend_position=:outertop)
plot!(title=L"He cooling: Q_{target} = %$(q/1e6) MW; $T_{inlet}$ = %$(T_inlet) K, $f_{circular}$ = %$(f_circular)")

plot!(titlefontsize=10, ylabelfontsize=15,xlabelfontsize=15, size=(500,500))


##%
plot()

δ = collect(0.001:0.001:0.01)
T_inlet = 273 + 200
α = collect(1:0.01:9.0).*1e4
α_,δ_ =  meshgrid(α,δ)
κ = 164
q = 2e7
f_circular = 1.07
T_max = @. (q/κ*δ_ + q/α_ + T_inlet) * f_circular

p3= contourf(δ,α, T_max;colorbar=true, clim=clims,aspect=1)
# cs = contour!(δ,α, T_max; levels=1200:1200, linewidth=3, contour_labels = true, label="1200K",fontsize=10, color = cgrad([:red,:red,:red]))
cs = contour!(δ,α, T_max; levels=1800:1800, linewidth=5, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="ITER baseline: 1800K")
# cs = contour!(δ,α, T_max; levels=1500:1500, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="1500K")
# cs = contour!(δ,α, T_max; levels=2000:2000, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="2000K")
plot!(;xlabel=L"\delta_{surface \rightarrow pipe} [m]")
plot!(;ylabel=L" heat transfert coeff. $[W.m^{-2}.K^{-1}]$")
#vline!([0.006];color=:green, linewidth=2,label="ITER monoblock design [Bonin]")
vline!([0.006];color=:blue,linestyle=:dash, linewidth=5,label="ITER monoblock design [Morshedy]")
hline!([0.026e6];linewidth=5,label="He: smooth Cu tube[Baxi]")
# hline!([0.04e6];linewidth=2,label="He: jet Cu tube [Baxi]")
hline!([0.05e6];linewidth=5,label="He: jet flow []")
hline!([0.08e6];color=:blue, linewidth=5,linestyle=:dot,label="Water: ITER")
plot!(;legend=nothing,legend_position=legend_position=:outertop)
plot!(title=L"He cooling: Q_{target} = %$(q/1e6) MW; $T_{inlet}$ = %$(T_inlet) K, $f_{circular}$ = %$(f_circular)")

plot!(titlefontsize=10, ylabelfontsize=15,xlabelfontsize=20, size=(600,600))

plot()

δ = collect(0.001:0.001:0.01)
T_inlet = 273 + 150
α = collect(1:0.01:9.0).*1e4
α_,δ_ =  meshgrid(α,δ)
κ = 164
q = 2e7
f_circular = 1.07
T_max = @. (q/κ*δ_ + q/α_ + T_inlet) * f_circular

p3= contourf(δ,α, T_max;colorbar=true, clim=clims,aspect=1)
# cs = contour!(δ,α, T_max; levels=1200:1200, linewidth=3, contour_labels = true, label="1200K",fontsize=10, color = cgrad([:red,:red,:red]))
cs = contour!(δ,α, T_max; levels=1800:1800, linewidth=5, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="ITER baseline: 1800K")
# cs = contour!(δ,α, T_max; levels=1500:1500, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="1500K")
# cs = contour!(δ,α, T_max; levels=2000:2000, linewidth=3, contour_labels = true, fontsize=10, color = cgrad([:blue,:blue,:blue]), label="2000K")
plot!(;xlabel=L"\delta_{surface \rightarrow pipe} [m]")
plot!(;ylabel=L" heat transfert coeff. $[W.m^{-2}.K^{-1}]$")
#vline!([0.006];color=:green, linewidth=2,label="ITER monoblock design [Bonin]")
vline!([0.006];color=:blue,linestyle=:dash, linewidth=5,label="ITER monoblock design [Morshedy]")
hline!([0.026e6];linewidth=5,label="He: smooth Cu tube[Baxi]")
# hline!([0.04e6];linewidth=2,label="He: jet Cu tube [Baxi]")
hline!([0.05e6];linewidth=5,label="He: jet flow []")
hline!([0.08e6];color=:blue, linewidth=5,linestyle=:dot,label="Water: ITER")
plot!(;legend=nothing,legend_position=legend_position=:outertop)
plot!(title=L"He cooling: Q_{target} = %$(q/1e6) MW; $T_{inlet}$ = %$(T_inlet) K, $f_{circular}$ = %$(f_circular)")

plot!(titlefontsize=10, ylabelfontsize=15,xlabelfontsize=20, size=(600,600))
