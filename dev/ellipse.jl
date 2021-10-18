using Plots
using Symbolics
using ForwardDiff 

# ellipse

a = 0.25 
b = 0.125 
xc = 0.0
yc = 1.0
function rotation_mat(x) 
    [cos(x) -sin(x); sin(x) cos(x)] 
end

α = 0.2 * π
θ = range(0.0, stop=2.0 * π, length=100)
xe(t) = a * cos(t) 
ye(t) = b * sin(t) 
ze(t, p) = rotation_mat(p[3]) * [xe(t); ye(t)] + p[1:2]
xt = xe.(θ)
yt = ye.(θ)
zt = hcat([ze(t, [xc; yc; α]) for t in θ]...)

plot(zt[1, :], zt[2, :], aspect_ratio=:equal, width=2.0, color=:orange)

# segment 
f(x) = 0.0 * x
t = range(-1.0, stop=1.0, length=100)

plot!(t, f.(t), width=2.0, color=:black)

d(w) = sqrt((ze(w[1], [xc; yc; α])[1] - w[2])^2.0 + (ze(w[1], [xc; yc; α])[2] - f(w[2]))^2.0)
d(w) = (ze(w[1], [xc; yc; α])[1] - w[2])^2.0 + (ze(w[1], [xc; yc; α])[2] - f(w[2]))^2.0

@variables w[1:2]

# dsym = d(w)
# dd = Symbolics.gradient(dsym, w)
# ddd = Symbolics.jacobian(dd, w)

# d(a, b) = d([a; b])
# plt = surface(-2.0:0.1:2.0, -2.0:0.1:2.0, d)
# plot!(plt[1], camera = (90,90))

# minimum distance 
function minimum_distance_point(w, p) 
    # distance
    d(w) = (ze(w[1], p)[1] - w[2])^2.0 + (ze(w[1], p)[2] - f(w[2]))^2.0
    r(w) = ForwardDiff.gradient(d, w) 
    rw(w) = ForwardDiff.hessian(d, w) 
    # initialization
    a = w
    # Newton
    for i = 1:10
        norm(r(a), 1) < 1.0e-8 && break
        # step
        Δ = -rw(a) \ r(a)
        α = 1.0
        a_cand = a + α * Δ
        # linesearch 
        iter = 0 
        while norm(r(a_cand), 1) >= norm(r(a), 1) && iter < 100
            α *= 0.5 
            a_cand = a + α * Δ
            iter += 1
        end
        a = a_cand
        @show norm(r(a), 1)
    end
    # check min-max
    pt1_surf = [a[2], f(a[2])]

    pt1_ellipse = ze(a[1], p)
    dis1 = norm(pt1_ellipse - pt1_surf) 

    # initialization
    a = w + [π; 0.0]
    # Newton
    for i = 1:10
        norm(r(a), 1) < 1.0e-8 && break
        # step
        Δ = -rw(a) \ r(a)
        α = 1.0
        a_cand = a + α * Δ
        # linesearch 
        iter = 0 
        while norm(r(a_cand), 1) >= norm(r(a), 1) && iter < 100
            α *= 0.5 
            a_cand = a + α * Δ
            iter += 1
        end
        a = a_cand
        @show norm(r(a), 1)
    end
    # check min-max
    pt2_surf = [a[2], f(a[2])]

    pt2_ellipse = ze(a[1], p)
    dis2 = norm(pt2_ellipse - pt2_surf) 

    if dis1 < dis2 
        dis = dis1 
        return pt1_ellipse, pt1_surf
    else
        dis = dis2 
        return pt2_ellipse, pt2_surf
    end
end

pt_ellipse, pt_surface = minimum_distance_point([0.0 * π; 0.0], [xc; yc; α])

plot!([pt_ellipse[1]; pt_surface[1]], [pt_ellipse[2]; pt_surface[2]], color=:red)

