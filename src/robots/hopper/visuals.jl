# Visualization
function visualize!(vis, model::Hopper, q; Δt = 0.1, fixed_camera=true)

    r_body = model.body_radius
    r_foot = model.foot_radius
    r_leg = 0.5 * r_foot

    default_background!(vis)

    setobject!(vis["body"], Sphere(Point3f0(0),
        convert(Float32, model.body_radius)),
        MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

    setobject!(vis["foot"], Sphere(Point3f0(0),
        convert(Float32, model.foot_radius)),
        MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, 1.0)))

    n_leg = 100
    for i = 1:n_leg
        setobject!(vis["leg$i"], Sphere(Point3f0(0),
            convert(Float32, r_leg)),
            MeshPhongMaterial(color = RGBA(0, 0, 0, 1.0)))
    end

    p_leg = [zeros(3) for i = 1:n_leg]
    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    for t = 1:length(q)
        p_body = [q[t][1], 0.0, q[t][2]]
        p_foot = [kinematics_foot(model, q[t])[1], 0.0, kinematics_foot(model, q[t])[2]]

        q_tmp = Array(copy(q[t]))
        r_range = range(0, stop = q[t][4], length = n_leg)
        for i = 1:n_leg
            q_tmp[4] = r_range[i]
            p_leg[i] = [kinematics_foot(model, q_tmp)[1], 0.0, kinematics_foot(model, q_tmp)[2]]
        end
        q_tmp[4] = q[t][4]
        p_foot = [kinematics_foot(model, q_tmp)[1], 0.0, kinematics_foot(model, q_tmp)[2]]

        z_shift = [0.0; 0.0; 0.0]

        MeshCat.atframe(anim, t) do
            settransform!(vis["body"], Translation(p_body + z_shift))
            settransform!(vis["foot"], Translation(p_foot + z_shift))

            for i = 1:n_leg
                settransform!(vis["leg$i"], Translation(p_leg[i] + z_shift))
            end
        end
    end

    if fixed_camera 
        settransform!(vis["/Cameras/default"],
        compose(Translation(0.0, -50.0, -1.0),LinearMap(RotZ(-pi / 2.0))))
        setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 25)
    end

    MeshCat.setanimation!(vis, anim)
end