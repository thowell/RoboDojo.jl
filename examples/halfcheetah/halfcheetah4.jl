# PREAMBLE

# PKG_SETUP

# ## Setup

using RoboDojo

# ## Initial conditions
q1 = nominal_configuration(halfcheetah4)
v1 = zeros(halfcheetah4.nq)
q1[2] += 0.25
q1[3] += 0.0 * π

# ## Time
h = 0.01
T = 1000

# ## policy 


# ## Simulator
s = Simulator(halfcheetah4, T, h=h)

# ## Simulate
simulate!(s, q1, v1)

# ## Visualizer
vis = Visualizer()
render(vis)

# ## Visualize
visualize!(vis, s)

# ## simulate policy
d = Simulator(halfcheetah4, 1, h=h)

x_hist = [[q1; q1 - v1 * h]]
u_hist = Vector{Float64}[]

for t = 1:T
    push!(u_hist, [0;0;0; 5.0e1 * randn(6)])
    y = zeros(2 * halfcheetah4.nq)
    dynamics(d, y, x_hist[end], u_hist[end], zeros(0))
    push!(x_hist, y)
end

s = Simulator(halfcheetah4, T, h=h)
for i = 1:T
    q = x_hist[i][1:halfcheetah4.nq]
    v = x_hist[i][halfcheetah4.nq .+ (1:halfcheetah4.nq)]
    set_state!(s, q, v, i)
end
visualize!(vis, s)


u_hist