function build_robot!(vis, model::Particle;
    i=1,
    r=0.1,
    tl=1.0,
    particle_color=Colors.RGBA(0.0, 0.0, 0.0, tl))

    setobject!(vis["particle_$i"],
        GeometryBasics.Sphere(Point3f0(0),
        convert(Float32, r)),
        MeshPhongMaterial(color=particle_color))
end

function set_robot!(vis, model::Particle, q; i=1)
    settransform!(vis["particle_$i"],
        compose(Translation(q...), LinearMap(RotY(0.0))))
end

function visualize!(vis, model::Particle, q;
    i=1,
    r=0.1,
    tl=1.0,
    particle_color=Colors.RGBA(0.0, 0.0, 0.0, tl),
    Δt=0.1,
    fixed_camera=true)

    shift = [0.0; 0.0; r]

    default_background!(vis)

    build_robot!(vis, model,
        i=i,
        r=r,
        tl=tl,
        particle_color=particle_color)

    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    T = length(q)
    for t = 1:T-1
        MeshCat.atframe(anim, t) do
            set_robot!(vis, model, q[t] + shift, i=i)
        end
    end

    if fixed_camera
        MeshCat.settransform!(vis["/Cameras/default"],
        MeshCat.compose(MeshCat.Translation(0.0, -50.0, -1.0), MeshCat.LinearMap(Rotations.RotZ(-pi / 2.0))))
        setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 25)
    end

    MeshCat.setanimation!(vis, anim)
end
