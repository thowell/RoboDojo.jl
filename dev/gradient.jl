# PREAMBLE

# PKG_SETUP

# ## Setup

using RoboDojo 
using Plots

q1 = nominal_configuration(box)
q1 = [0.0; 0.5; 0.0]
v1 = zeros(box.nq)
u1 = zeros(box.nu) 

T = 1
h = 0.1

u = range(0.0, stop=2.0, length=1000)
κ = [10.0^(-i) for i = 2:8]

x = []
∂x∂u = [] 
κtol = []

for (i, κi) in enumerate(κ)
    xi = [] 
    ∂xi∂u = []
    κitol = []
    for u in ux
        s = Simulator(box, T, h=h,
            solver_opts=InteriorPointOptions(
                undercut=2.0,
                r_tol=1e-10,
                κ_tol=κi,  
                diff_sol=true,
                verbose=false))

        step!(s, q1, v1, [u; 0.0; 0.0], 1)

        push!(xi, s.traj.q[end][1]) 
        push!(∂xi∂u, s.grad.∂q3∂u1[end][1, 1])
        push!(κitol, s.ip.κ[1])
    end
    push!(x, xi) 
    push!(∂x∂u, ∂xi∂u)
    push!(κtol, κitol)
end

plt_x = plot(; xlabel="F", ylabel="Δx", label=""); 
for (i, xi) in enumerate(x) 
    plt_x = plot!(u, xi, label="κ=$(scn(κ[i]))", legend=:topleft) 
end
display(plt_x)

plt_∂x∂u = plot(; xlabel="F", ylabel="∂x∂F", label=""); 
for (i, ∂xi∂u) in enumerate(∂x∂u) 
    plt_∂x∂u = plot!(u, ∂xi∂u, label="κ=$(scn(κ[i]))", legend=:topleft) 
end
display(plt_∂x∂u)


# ## Visualizer
vis = Visualizer()
render(vis)

# ## Visualize
visualize!(vis, s)