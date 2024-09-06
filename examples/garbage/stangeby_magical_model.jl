#=
Author: Jerome Guterl (guterlj@fusion.gat.com)
 Company: General Atomics
 stangeby_magical_model.jl (c) 2024=#

# Some magic happens...

Tₑₜ = 10 .^ LinRange(-1.0, 2.0, 1000)

plot(Tₑ, fcool.(Tₑ);yscale=:log10, xscale=:log10, ylim=[1e-2,1.0])

Tsep = [4.0, 2.0, 2.5, 4.0, 2.0]
frad_outer = [0.83, 0.81, 0.80, 0.86, 0.96]
for i in 1:5 
scatter!([Tsep[i]], [1 .- frad_outer[i]], label = "$i")
end
plot!()
q_para = 1e8
L_para = 100.0
f_cond = 1.0
κ_0 = 2930.0
Rₜ = 1.0
Rᵤ = 1.0
α = 10 .^ LinRange(5,10,1000)
Tᵤ = 250.0
meshgrid(x, y) = first.(Iterators.product(x, y)), last.(Iterators.product(x, y))
using NLsolve
function get_Te_target(Tᵤ,Rₜ,Rᵤ,q_para,L_para,f_cond,κ_0; Te_init=10.0 )
G(Tₑₜ) = @. Tᵤ^(7/2) - (Tₑₜ[1]^(7 / 2) + 7 / 2 * Rₜ / Rᵤ * q_para * L_para * f_cond / κ_0 / (1.0 - fcool(Tₑₜ[1])))
    sol = nlsolve(G, [Te_init])
sol.zero[1]
end

a = [-1.171, 2.101, -3.212, 3.567, -1.477]
_fcool(T) = 1.0 - 10^(sum([a_ * (log10(T))^(i - 1) for (i, a_) in enumerate(a)]))
fcool(T) = T > 10.0 ? 1.0 - _fcool(10.0) : 1.0 - _fcool(T)

function get_frad(Tᵤ, Rₜ, Rᵤ, q_para, L_para, f_cond, κ_0; Te_init=10.0)
    Te_target = get_Te_target(Tᵤ, Rₜ, Rᵤ, q_para, L_para, f_cond, κ_0; Te_init=10.0)
    fcool(Te_target)
end

compute_Tₑᵤ(Tₑₜ, α) = @. (Tₑₜ^(7 / 2) + α / (1 - fcool(Tₑₜ)))^(2 / 7)
(Tₑₜ_,α_) =  meshgrid(Tₑₜ, α)
Tᵤ = compute_Tₑᵤ.(meshgrid(Tₑₜ, α)...)
Tᵤ  = interpolate((Tₑₜ, α), compute_Tₑᵤ.(meshgrid(Tₑₜ, α)...), (Gridded(Linear()), Gridded(Linear())))

α = 10 .^ LinRange(5, 10,6)
Tᵤ(Tₑₜ, α)
using Interpolations
plot()
Tₑₜ = 10 .^ LinRange(-1.0, 1.0, 100)

for α_ in 10 .^ LinRange(5, 10,6)
    plot!(1 .- fcool.(Tₑₜ), Tₑᵤ(Tₑₜ, α_ .+ 0 .* Tₑₜ); label="α = $α_")
end
plot!()


