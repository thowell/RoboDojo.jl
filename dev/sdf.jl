using Symbolics 
using ForwardDiff
using LinearAlgebra
using Plots
using ForwardDiff

struct Segment{T}
    domain::Tuple{T,T}
    func
    f 
    df
    ddf
    r 
    rx 
    rp 
    δp
    dir_norm
    δ_dir_norm
end

function Segment(domain, func) 
    @variables x p[1:2] 

    # surface
    f = func(x) 
    df = Symbolics.derivative(f, x) 
    ddf = Symbolics.derivative(df, x) 

    # distance (squared) 
    r = x - p[1] + (f - p[2]) * df
    rx = 1.0 + df * df + (f - p[2]) * ddf
    rp = [-1.0 -df]
    δp = [-rp ./ rx; -df * rp ./ rx]

    # fast methods
    f_func = eval(Symbolics.build_function(f, x))
    df_func = eval(Symbolics.build_function(df, x))
    ddf_func = eval(Symbolics.build_function(ddf, x))

    r_func = eval(Symbolics.build_function(r, x, p))
    rx_func = eval(Symbolics.build_function(rx, x, p))
    fp_func = eval(Symbolics.build_function(rp, x, p)[1])
    δp_func = eval(Symbolics.build_function(δp, x, p)[1])

    # surface norm 
    f2 = func(p[1]) 
    dir = [-Symbolics.derivative(f2, p[1]); 1.0]
    dir_norm = dir ./ norm(dir)
    dir_norm_func = eval(Symbolics.build_function(dir_norm, p)[1])
    δ_dir_norm = Symbolics.jacobian(dir_norm, p) 
    δ_dir_norm_func = eval(Symbolics.build_function(δ_dir_norm, p)[1])

    Segment(domain, func, 
        f_func, df_func, ddf_func, 
        r_func, rx_func, fp_func, 
        δp_func,
        dir_norm_func, δ_dir_norm_func)
end
# seg1 = Segment((-10.0, 10.0), x -> 0.5 * x^2.0)

seg1 = Segment((-2.0, 0.5), x -> 0.0 * x)
seg2 = Segment((0.5, 1.5), x -> x - 0.5)
# seg3 = Segment((1.5, 2.5), x -> 1.0)
# seg4 = Segment((2.5, 2.51), x -> 100.0 * x - 249.0)
# seg5 = Segment((2.51, 3.51), x -> 2.0)
# seg6 = Segment((3.51, 3.52), x -> -100.0 * x + 353.0)
# seg7 = Segment((3.52, 4.52), x -> 1.0)

surf = [seg1, seg2]#, seg3]#, seg4]#, seg5, seg6, seg7]

get_point(s::Segment, a) = [a; s.f(a)] 

p = [0.0; 1.0]

function minimum_distance_point(s::Segment, p) 
    # initialization
    if p[1] < s.domain[1] 
        a = s.domain[1] 
    elseif p[1] > s.domain[2] 
        a = s.domain[2] 
    else 
        a = p[1] 
    end
   
    # Newton
    for i = 1:100
        abs(s.r(a, p)) < 1.0e-8 && break
        # step
        Δ = -s.rx(a, p) \ s.r(a, p)
        α = 1.0
        a_cand = a + α * Δ
        # linesearch 
        iter = 0 
        while abs(s.r(a_cand, p)) >= abs(s.r(a, p)) && iter < 100
            α *= 0.5 
            a_cand = a + α * Δ
            iter += 1
        end
        a = a_cand
    end
        
    # solution
    if a < s.domain[1] 
        a = s.domain[1]
    elseif a > s.domain[2] 
        a = s.domain[2]
    end 

    # point and sensitivity
    return get_point(s, a), s.δp(a, p)
end

get_surf_dir(s::Segment, a) = [-s.df(a[1]); 1.0]

function surface_normal(s::Segment, p) 
    dir = get_surf_dir(s, p)
    dir_norm = dir ./ norm(dir)
    # dir_norm = s.dir_norm(p)

    θ = -0.5 * π
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    tan_norm = R * dir_norm

    return dir_norm, tan_norm
end

function rotation(ns)
	# world-frame normal
	nw = [0.0, 1.0]

	angle = atan(nw[2], nw[1]) - atan(ns[2], ns[1])

	# Rw->s
	[cos(angle) -sin(angle); sin(angle) cos(angle)]
end

function check_domain(surf::Vector{Segment{T}}, p; r=0.0) where T
    status = Bool[] 
    for s in surf 
        if p[1] >= s.domain[1] - r && p[1] <= s.domain[2] + r
            push!(status, true) 
        else 
            push!(status, false) 
        end
    end
    return status 
end

function contact_point(s::Segment, p) 
    pmin, δp = minimum_distance_point(s, p)
    d = norm(p - pmin)
    dir_norm, tan_norm = surface_normal(s, pmin)
    rot = rotation(dir_norm)
    δ_dir_norm = s.δ_dir_norm(pmin)
    return pmin, d, dir_norm, tan_norm, rot, δp, δ_dir_norm
end

function contact_point(surf::Vector{Segment{T}}, p; r=0.0) where T 
    domain = check_domain(surf, p, r=r)
    p_eval = []
    for (i, status) in enumerate(domain) 
        if status
            push!(p_eval, contact_point(surf[i], p))
        end
    end 
    min_idx = findmin([pt[2] for pt in p_eval])[2]
    return p_eval[min_idx]
end

function plot_segment(s::Segment; N=100) 
    x = range(s.domain[1], stop=s.domain[2], length=N)
    y = s.f.(x) 
    plot!(x, y, width=2.0, color=:black, label="") 
end

function visualize_segment!(vis, s::Segment; 
    idx = 1, color=RGBA(0.0, 0.0, 0.0, 1.0), linewidth=50, num_points=100)

    line_mat = LineBasicMaterial(color=color, linewidth=linewidth)
    domain = range(s.domain[1], stop=s.domain[2], length=num_points)
    points = [Point(x, 0.0, s.func(x)) for x in domain]
    setobject!(vis["traj_$idx"], MeshCat.Line(points, line_mat))
end

function visualize!(vis, surf::Vector{Segment{T}}; 
    colors=[RGBA(0.0, 0.0, 0.0, 1.0) for i = 1:length(surf)], linewidth=50, num_points=100) where T

    for (i, s) in enumerate(surf)
        visualize_segment!(vis, s, idx = i, color=colors[i], linewidth=linewidth, num_points=num_points)
    end
end

plt = plot(xlims=(-0.5, 4.0), ylims=(-1.0, 3.0), aspect_ratio=:equal)
for s in surf 
    plt = plot_segment(s)
end

# p = [2.0; 1.5]
# p = [0.0, 1.0]
# p = [0.0, 1.4]
p = [0.5, 1.0]
p = [1.5, 1.25]

plt = scatter!([p[1]], [p[2]], color=:red, label="")

p_eval = [contact_point(s, p) for s in surf]
for x in p_eval
    plt = plot!([p[1]; x[1][1]], [p[2]; x[1][2]], width=1.5, color=:cyan, label="")
end
display(plt)

for (i, pc) in enumerate(p_eval)
    @show x = pc[1][1] 
    @show y = pc[1][2]

    dir_norm = pc[3] 
    tan_norm = pc[4]

    scale = 0.25
    plt = plot!([x; x + scale * dir_norm[1]], [y; y + scale * dir_norm[2]], width=1.5, color=:magenta, label="")
    plt = plot!([x; x + scale * tan_norm[1]], [y; y + scale * tan_norm[2]], width=1.5, color=:orange, label="")

end

display(plt)

function generate_polynomial(seg1::Segment{T}, seg2::Segment{T}; ϵ=0.01) where {T}
    # segments
    seg1.domain[2]
    seg2.domain[1]

    # offsets
    x⁺ = seg1.domain[2] + ϵ
    x⁻ = seg1.domain[2] - ϵ

    # end points and their derivatives
    p1x = seg1.func(x⁻) 
    dp1x = ForwardDiff.derivative(seg1.func, x⁻)
    p2x = seg2.func(x⁺) 
    dp2x = ForwardDiff.derivative(seg2.func, x⁺)

    # fit polynomial parameters
    b = [p1x; dp1x; p2x; dp2x]
    A = [1.0 x⁻ x⁻^2.0 x⁻^3.0; 
        0.0 1.0 2.0 * x⁻ 3.0 * x⁻^2.0; 
        1.0 x⁺ x⁺^2.0 x⁺^3.0;
        0.0 1.0 2.0 * x⁺ 3.0 * x⁺^2.0]
    param = A \ b

    # polynomial
    function poly(y) 
        param[1] + param[2] * y + param[3] * y^2.0 + param[4] * y^3.0
    end
    @assert abs(poly(x⁻) - seg1.func(x⁻)) < 1.0e-8
    @assert abs(poly(x⁺) - seg2.func(x⁺)) < 1.0e-8

    # polynomial segment
    Segment((x⁻, x⁺), poly)
end

# seg_poly = generate_polynomial(seg1, seg2)

function polynomial_blend(seg1::Segment{T}, seg2::Segment{T}; ϵ=0.01) where {T}
    # polynominal segment
    seg_poly = generate_polynomial(seg1, seg2, ϵ=ϵ)

    # update segments
    seg1_update = Segment((seg1.domain[1], seg_poly.domain[1]), seg1.func)
    seg2_update = Segment((seg_poly.domain[2], seg2.domain[2]), seg2.func)

    return seg1_update, seg_poly, seg2_update 
end

seg1_update, seg_poly, seg2_update = polynomial_blend(seg1, seg2)

function polynomial_blend(surf::Vector{Segment{T}}; ϵ=0.01) where {T}
    @assert length(surf) > 1
    seg1_update, seg_poly, seg2_update = polynomial_blend(surf[1], surf[2], ϵ=ϵ)
    surf_update = [seg1_update, seg_poly, seg2_update]
    if length(surf) > 2 
        for i = 2:length(surf)-1
            seg1_update, seg_poly, seg2_update = polynomial_blend(surf_update[end], surf[i+1], ϵ=ϵ)
            pop!(surf_update) 
            push!(surf_update, seg1_update) 
            push!(surf_update, seg_poly) 
            push!(surf_update, seg2_update) 
        end 
    end
    return surf_update
end

surf_update = polynomial_blend(surf, ϵ=0.1)

# plt = plot(xlims=(1.45, 1.55), ylims=(0.95, 1.05), aspect_ratio=:equal)
plt = plot(xlims=(-0.5, 4.0), ylims=(-1.0, 2.5), aspect_ratio=:equal)
# plt = plot(xlims=(3.45, 3.55), ylims=(0.95, 1.05), aspect_ratio=:equal)
plt = plot(xlims=(0.0, 1.0), ylims=(-0.1, 0.25), aspect_ratio=:equal)

for s in surf_update
    plt = plot_segment(s)
end
display(plt)

p = [1.4; 1.25]
plt = scatter!([p[1]], [p[2]], color=:red, label="")

p_eval = [contact_point(s, p) for s in surf_update]
for x in p_eval
    plt = plot!([p[1]; x[1][1]], [p[2]; x[1][2]], width=1.5, color=:cyan, label="")
end
display(plt)

for (i, pc) in enumerate(p_eval)
    @show x = pc[1][1] 
    @show y = pc[1][2]

    dir_norm = pc[3] 
    tan_norm = pc[4]

    scale = 0.25
    plt = plot!([x; x + scale * dir_norm[1]], [y; y + scale * dir_norm[2]], width=2.0, color=:magenta, label="")
    plt = plot!([x; x + scale * tan_norm[1]], [y; y + scale * tan_norm[2]], width=2.0, color=:orange, label="")
end
display(plt)

# # test rotation
# p_eval[3][1]
# rotation(p_eval[6][3])' * [1.0; 0.0] # map vector in segment to world frame