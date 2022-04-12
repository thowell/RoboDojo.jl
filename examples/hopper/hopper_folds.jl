# PREAMBLE

# PKG_SETUP

# ## Setup

using RoboDojo 
using BenchmarkTools 
using Folds

# ## Initial conditions
q1 = nominal_configuration(hopper)
v1 = zeros(hopper.nq)

# ## Time 
T = 100
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

# ## Simulate

using Symbolics
n = 4
@variables a[1:n] b[1:n]
function cone_product(u, v)
    [transpose(u) * v; u[1] * v[2:end] + v[1] * u[2:end]]
end

using LinearAlgebra

function jacobian(u) 
    n = length(u)
    U = u[1] * Array(Diagonal(ones(n)))
    U[2:end,1] = u[2:end]
    U[1,2:end] = u[2:end]
    return U
end

function cp_inverse(u)
    n = length(u)
    α = -1/u[1]^2 * norm(u[2:end])^2
    β = 1 / (1 + α)
    # S1 = zeros(n,n)
    # S1[end,1:end-1] = u[end:-1:2]/u[1]

    S1 = [
            zeros(n-1, n); 
            u[end:-1:2]'/u[1] 0.0;
    ]
    # S2 = zeros(n,n)
    # S2[1:end-1,end] = u[end:-1:2]/u[1]

    S2 = [
        zeros(n-1, n-1)  u[end:-1:2]/u[1];
        zeros(1, n)
    ]
    P = zeros(n,n)
    for i = 1:n
        P[end-i+1,i] = 1
    end
    U = jacobian(u)
    V = P * 1/u[1] * U * P
    Vi = (I - S1) * (I - β * (S2 * (I - S1)))
    Ui = P * u[1] * I(n) * Vi * P
    return Ui
end

c = cone_product(a, b)
ca = Symbolics.jacobian(c, a)
cb = Symbolics.jacobian(c, b)
cai = cp_inverse(b)

jacobian(ones(n)) * cp_inverse(ones(n))
ca * cai
cp_inverse(ones(n))

n = 3
x = rand(n)
u = [1;rand(n-1)]#[1e-10,1,0]
U = jacobian(u)
U * inverse(u) * x - x