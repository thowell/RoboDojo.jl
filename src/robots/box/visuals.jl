function build_robot!(vis, model::Box;
    i=1,
    r=0.1,
    tl=1.0,
    box_color=Colors.RGBA(0.0, 0.0, 0.0, tl))

    setobject!(vis["box_$i"], GeometryBasics.Rect(
        Vec(-1.0 * r, -0.005, -1.0 * r),
        Vec(2.0 * r, 0.01, 2.0 * r)),
        MeshPhongMaterial(color=box_color))
end

function set_robot!(vis, model::Box, q; i=1)
    settransform!(vis["box_$i"],
        compose(Translation(q[1], 0.001 * i, q[2]), LinearMap(RotY(-1.0 * q[3]))))
end

function visualize!(vis, model::Box, q;
    i=1,
    r=0.1,
    tl=1.0,
    box_color=Colors.RGBA(0.0, 0.0, 0.0, tl),
    Δt=0.1,
    fixed_camera=true)

    default_background!(vis)

    build_robot!(vis, model,
        i=i,
        r=r,
        tl=tl,
        box_color= box_color) 

    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    T = length(q)
    for t = 1:T-1
        MeshCat.atframe(anim, t) do
            set_robot!(vis, model, q[t], i=i)
        end
    end

    if fixed_camera
        MeshCat.settransform!(vis["/Cameras/default"],
        MeshCat.compose(MeshCat.Translation(0.0, -50.0, -1.0), MeshCat.LinearMap(Rotations.RotZ(-pi / 2.0))))
        setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 25)
    end

    MeshCat.setanimation!(vis, anim)
end