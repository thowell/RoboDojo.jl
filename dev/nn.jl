using Plots, Symbolics

# environment points
pa = [0.0; 0.0] 

pb = [1.0; 0.0] 

pc = [1.0; 1.0] 

pd = [2.0; 1.0] 

pe = [2.0; 0.0] 

pf = [3.0; 0.0] 

pts_env = [pa, pb, pc, pd, pe, pf]
plot(hcat(pts_env...)[1, :], hcat(pts_env...)[2, :], aspect_ratio=:equal, color=:magenta, width=2.0)

function points_between(p1, p2; N = 100)  
    dir = p2 - p1 
    pts = [p1 + t * dir for t in range(0.0, stop=1.0, length=N)]
    return pts 
end

N = 2
pts_ab = points_between(pa, pb, N=N)
pts_bc = points_between(pb, pc, N=N)
pts_cd = points_between(pc, pd, N=N)
pts_de = points_between(pd, pe, N=N)
pts_ef = points_between(pe, pf, N=N)

pts_data = [pts_ab..., pts_bc...]#, pts_cd..., pts_de..., pts_ef...]

plot!(hcat(pts_data...)[1, :], hcat(pts_data...)[2, :], aspect_ratio=:equal, color=:orange, width = 1.0)

# data set
X = hcat(pts_data...)[1, :]
Y = hcat(pts_data...)[2, :]

# model
l_input = 1 
l_1 = 2 
l_2 = 2 
l_output = 1
p = l_1 * l_input + l_1 + l_2 * l_1 + l_2 + l_output * l_2 + l_output 

function model(x, θ) 
    W1 = reshape(θ[1:(l_1 * l_input)], l_1, l_input)
    b1 = θ[(l_1 * l_input) .+ (1:l_1)]

    W2 = reshape(θ[(l_1 * l_input + l_1) .+ (1:(l_2 * l_1))], l_2, l_1)
    b2 = θ[(l_1 * l_input + l_1 + l_2 * l_1) .+ (1:l_2)] 

    W3 = reshape(θ[(l_1 * l_input + l_1 + l_2 * l_1 + l_2) .+ (1:(l_output * l_2))], l_output, l_2)
    b3 = θ[(l_1 * l_input + l_1 + l_2 * l_1 + l_2 + l_output * l_2) .+ (1:l_output)] 

    z1 = W1 * x + b1 
    return [sum(z1)]
    o1 = max.(z1, 0.0) 
     
    z2 = W2 * o1 + b2 
    o2 = max.(z2, 0.0) 

    z3 = W3 * o2 + b3 

    return z3 
end

@variables x[1:1], θ[1:p]

m = model(x, θ)
dm = vec(Symbolics.jacobian(m, θ))

m_func = eval(Symbolics.build_function(m, x, θ)[1])
dm_func = eval(Symbolics.build_function(dm, x, θ)[1]) 

m_func([2.0], θ0)
function eval_obj(X, Y, θ) 
    N = length(X) 
    J = 0.0 
    for i = 1:N
        J += (m_func(X[i:i], θ) - Y[i:i])[1]^2.0 
    end 
    return 0.5 * J / N
end

function eval_grad(X, Y, θ) 
    N = length(X) 
    ∇J = zero(θ)
    for i = 1:N
        ∇J .+= dm_func(X[i:i], θ) * (m_func(X[i:i], θ) - Y[i:i])[1]
    end 
    return ∇J ./ N
end

θ0 = 0.01 * randn(p)
eval_obj(X, Y, θ0)
eval_grad(X, Y, θ0)


function optimize!(X, Y, θ; max_iter=1000, α=1.0e-6, epoch=10) 
    println("obj (iter 1): $(eval_obj(X, Y, θ))")

    for i = 1:max_iter 
        ∇J = eval_grad(X, Y, θ) 
        θ .-= α * ∇J 

        i % epoch == 0 && (println("obj (iter $i): $(eval_obj(X, Y, θ))"))
    end
    return θ 
end

θ = 0.01 * randn(p)
θ = optimize!(X, Y, θ, α=1.0e-3, max_iter=100000, epoch=1000);

X_pred = range(0, stop=3.0, length = 100) 
Y_pred = [model(X_pred[i:i], θ)[1] for i = 1:length(X_pred)] 

plot(X_pred, Y_pred, aspect_ratio=:equal, color=:cyan, width = 1.0)

model([1.0], θ0)
