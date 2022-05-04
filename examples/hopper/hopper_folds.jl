# PREAMBLE

# PKG_SETUP

# ## Setup

using RoboDojo 
using BenchmarkTools 
using Folds

@show Threads.nthreads()

# ## Initial conditions
q1 = nominal_configuration(hopper)
v1 = zeros(hopper.nq)

# ## Time 
T = 1
h = 0.01

# ## Simulator 
s = Simulator(hopper, T, h=h)

# ## Simulate
@benchmark simulate!(s, q1, v1) setup=(s=s, q1=q1, v1=v1)

# ## Visualizer
vis = Visualizer()
render(vis)

# ## Visualize
visualize!(vis, s)

# ## for loop 
function sim_forloop!(S, Q, V, N) 
    for i = 1:N 
        simulate!(S[i], Q[i], V[i])
    end 
    return 
end

function sim_threads!(S, Q, V, N) 
    Threads.@threads for i = 1:N 
        simulate!(S[i], Q[i], V[i])
    end 
    return 
end

N = 16
S = [deepcopy(s) for i = 1:N]
Q = [q1 for i = 1:N] 
V = [v1 for i = 1:N]
X = [(S[i], Q[i], V[i]) for i = 1:N]

sim_forloop!(S, Q, V, N)
@benchmark sim_forloop!(S, Q, V, N) setup=(S=S, Q=Q, V=V, N=N)

sim_threads!(S, Q, V, N)
@benchmark sim_threads!(S, Q, V, N) setup=(S=S, Q=Q, V=V, N=N)
visualize!(vis, S[1])

Folds.map(simulate!, S, Q, V, Folds.ThreadedEx())
@benchmark Folds.map(simulate!, S, Q, V) setup=(simulate! =simulate!, S=S, Q=Q, V=V)
@benchmark Folds.map(simulate!, S, Q, V, Folds.ThreadedEx()) setup=(simulate! =simulate!, S=S, Q=Q, V=V)

