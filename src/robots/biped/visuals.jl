function build_robot!(vis::Visualizer, model::Biped; 
    r=0.035, r_contact=r * 8.0 / 7.0, color_opacity=1.0)

    r = convert(Float32, r) 
    r_contact = convert(Float32, r_contact)
    
    body_mat = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, color_opacity))
    contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, color_opacity))

    default_background!(vis)

    torso = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_torso), r)
    setobject!(vis[:robot]["torso"], torso, body_mat)

    thigh_1 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_thigh1), r)
    setobject!(vis[:robot]["thigh1"], thigh_1, body_mat)

    calf_1 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_calf1), r)
    setobject!(vis[:robot]["calf1"], calf_1, body_mat)

    foot_1 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.lt_foot1 + model.lh_foot1), r)
    setobject!(vis[:robot]["foot1"], foot_1, body_mat)

    thigh_2 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_thigh2), r)
    setobject!(vis[:robot]["thigh2"], thigh_2, body_mat)

    calf_2 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_calf2), r)
    setobject!(vis[:robot]["calf2"], calf_2, body_mat)

    foot_2 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.lt_foot2 + model.lh_foot2), r)
    setobject!(vis[:robot]["foot2"], foot_2, body_mat)

    setobject!(vis[:robot]["heel1"], Sphere(Point3f0(0.0), r_contact), contact_mat)
    setobject!(vis[:robot]["heel2"], Sphere(Point3f0(0.0), r_contact), contact_mat)
    setobject!(vis[:robot]["toe1"], Sphere(Point3f0(0.0), r_contact), contact_mat)
    setobject!(vis[:robot]["toe2"], Sphere(Point3f0(0.0), r_contact), contact_mat)
    setobject!(vis[:robot]["knee1"], Sphere(Point3f0(0.0), r), body_mat)
    setobject!(vis[:robot]["knee2"], Sphere(Point3f0(0.0), r), body_mat)
    setobject!(vis[:robot]["hip"], Sphere(Point3f0(0.0), r), body_mat)
    setobject!(vis[:robot]["torso_top"], Sphere(Point3f0(0.0), r), body_mat)

    return true
end

function set_robot!(vis::Visualizer, model::Biped, q;
   r=0.035, r_contact=r * 8.0 / 7.0)

    p_shift = [0.0; 0.0; r_contact]

    k_hip = kinematics_body(model, q, mode=:hip)
    p_hip = [k_hip[1]; 0.0; k_hip[2]] + p_shift

    k_torso = kinematics_body(model, q, mode=:ee)
    p_torso = [k_torso[1], 0.0, k_torso[2]] + p_shift

    k_thigh_1 = kinematics_thigh(model, q, leg=:leg1, mode=:ee)
    p_thigh_1 = [k_thigh_1[1], 0.0, k_thigh_1[2]] + p_shift

    k_calf_1 = kinematics_calf(model, q, leg=:leg1, mode=:ee)
    p_calf_1 = [k_calf_1[1], 0.0, k_calf_1[2]] + p_shift

    k_toe_1 = kinematics_foot(model, q, leg=:leg1, mode=:toe)
    p_toe_1 = [k_toe_1[1], 0.0, k_toe_1[2]] + p_shift

    k_heel_1 = kinematics_foot(model, q, leg=:leg1, mode=:heel)
    p_heel_1 = [k_heel_1[1], 0.0, k_heel_1[2]] + p_shift

    k_thigh_2 = kinematics_thigh(model, q, leg=:leg2, mode=:ee)
    p_thigh_2 = [k_thigh_2[1], 0.0, k_thigh_2[2]] + p_shift

    k_calf_2 = kinematics_calf(model, q, leg=:leg2, mode=:ee)
    p_calf_2 = [k_calf_2[1], 0.0, k_calf_2[2]] + p_shift

    k_toe_2 = kinematics_foot(model, q, leg=:leg2, mode=:toe)
    p_toe_2 = [k_toe_2[1], 0.0, k_toe_2[2]] + p_shift

    k_heel_2 = kinematics_foot(model, q, leg=:leg2, mode=:heel)
    p_heel_2 = [k_heel_2[1], 0.0, k_heel_2[2]] + p_shift

    settransform!(vis[:robot]["thigh1"], cable_transform(p_hip, p_thigh_1))
    settransform!(vis[:robot]["calf1"], cable_transform(p_thigh_1, p_calf_1))
    settransform!(vis[:robot]["foot1"], cable_transform(p_toe_1, p_heel_1))

    settransform!(vis[:robot]["thigh2"], cable_transform(p_hip, p_thigh_2))
    settransform!(vis[:robot]["calf2"], cable_transform(p_thigh_2, p_calf_2))
    settransform!(vis[:robot]["foot2"], cable_transform(p_toe_2, p_heel_2))

    settransform!(vis[:robot]["torso"], cable_transform(p_hip, p_torso))
    settransform!(vis[:robot]["heel1"], Translation(p_heel_1))
    settransform!(vis[:robot]["heel2"], Translation(p_heel_2))
    settransform!(vis[:robot]["toe1"], Translation(p_toe_1))
    settransform!(vis[:robot]["toe2"], Translation(p_toe_2))
    settransform!(vis[:robot]["knee1"], Translation(p_thigh_1))
    settransform!(vis[:robot]["knee2"], Translation(p_thigh_2))
    settransform!(vis[:robot]["hip"], Translation(p_hip))
    settransform!(vis[:robot]["torso_top"], Translation(p_torso))

    return true
end

function visualize!(vis, model::Biped, q; 
	Δt=0.1, r=0.04, r_contact=0.04, color_opacity=1.0, fixed_camera=true)

	build_robot!(vis, model, r=r, r_contact=r_contact, color_opacity=color_opacity)

	anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

	for (t, qt) in enumerate(q)
		MeshCat.atframe(anim, t) do
			set_robot!(vis, model, qt)
		end
	end

	MeshCat.setanimation!(vis, anim)

    if fixed_camera
        settransform!(vis["/Cameras/default"],
            compose(Translation(0.0, -50.0, -1.0),LinearMap(RotZ(-pi / 2.0))))
        setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 25)
    end

    return true
end