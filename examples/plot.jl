plot!(v, fraction; label="Ar")
v = []
nu = 1e20
qtarget = 1e7
Tup = 200.0
fraction = collect(0.001:0.001:0.2)
for f in fraction
    push!(v, sqrt((nu * Tup * V_legyel__ADAS(0.0, Tup, f, :Ne))^2 + qtarget^2) * 1e-6 * 1.35e-3 * (2 * pi))
end
plot!(v, fraction; label="Ne")
using Plots
gr()
#%%
using BoundaryPlasmaModels
using Colors
colors = palette(:tab10)
l = @layout [a b; b c]
p = plot(; layout=l)
imps = [:Ne, :Ar, :Kr, :Xe]
nu = [5e19, 1e20]
Te = [100.0, 200.0]
fraction = 10 .^ LinRange(-4:0.05:-1)
for (i, imp) in enumerate(imps)
    for (j, (nu_, Tup_)) in enumerate(Base.Iterators.product(nu, Te))
        q_rad = []
        for f in fraction
            push!(q_rad, nu_ * Tup_ * BoundaryPlasmaModels.V_legyel_ADAS(0.0, Tup_, f, imp) / 1e9)
        end
        plot!(fraction, q_rad; label="$nu_ , $Tup_", linewidth=3.0, color=colors[j], subplot=i, title="$imp", ylim=[0, 6], xlim=[0, 0.1])
    end

end
p

path = "/Users/jeromeguterl/Dropbox/Apps/Overleaf/Divertor design and analysis/figures/"
Plots.PGFPlotsX.export2tikz(path * "lengyel_integrale.tex",p;relpath="figures") #be aware of not override

#%%


plot!([1], [0], label="Lengyel model: GASC", color="black")
plot!([1], [0], linestyle=:dash, label="Lengyel model: 1D/2D informed", color="black")
plot!(xlabel=L"$P_{SOL} B_θ / R_{omp}#$",
    ylabel=L"$Z_{eff}#$ at separatrix upstream",
    legend_position=:outertop,
    frame=true, xlim=[1, 25],
    ylim=[1.0, 5.0],
    aspect_ratio=(24) / (4),
    xlabelfontsize=20,
    ylabelfontsize=20,
    size=(600, 600),
    legendfontsize=12,
    legend_columns=1,
    xtickfontsize=15,
    ytickfontsize=15
)
pgfplotsx()

p = plot()
x = collect(0:0.001:1.0)
y = @. 1 - sqrt(1 - x^2)
plot!(y, x, ylabel="Δq/qᵤ", xlabel="q_rad/qᵤ")
path = "/Users/jeromeguterl/Dropbox/Apps/Overleaf/Divertor design and analysis/figures/"
#Plots.PGFPlotsX.export2tikz(path * "delta_q.tex",p;relpath="figures") #be aware of not override


