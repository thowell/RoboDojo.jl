@testset "Robots: hopper" begin
    q0 = nominal_configuration(hopper) 

    @test norm(RobotDojo.kinematics_body(hopper, q0) - q0[1:2]) < 1.0e-8 
    @test norm(RobotDojo.kinematics_body_jacobian(hopper, q0) - [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0]) < 1.0e-8
    @test norm(RobotDojo.kinematics_foot(hopper, q0) - [0.0; 0.5]) < 1.0e-8
    @test norm(RobotDojo.kinematics_foot_jacobian(hopper, q0) - [1.0 0.0 (q0[4] * cos(q0[3])) sin(q0[3]); 0.0 1.0 (q0[4] * sin(q0[3])) (-1.0 * cos(q0[3]))]) < 1.0e-8
    @test lagrangian(hopper, zeros(hopper.nq), zeros(hopper.nq)) ≈ 0.0
    @test norm(RobotDojo.input_jacobian(hopper, q0) - transpose([0.0 -sin(q0[3]); 0.0 cos(q0[3]); 1.0 0.0; 0.0 1.0])) < 1.0e-8
    @test norm(RobotDojo.contact_jacobian(hopper, q0) - [1.0 0.0 0.0 0.0; # note that body torque due to friction is added in residual
        0.0 1.0 0.0 0.0;
        1.0 0.0 (q0[4] * cos(q0[3])) sin(q0[3]);
        0.0 1.0 (q0[4] * sin(q0[3])) (-1.0 * cos(q0[3]));
        0.0 0.0 0.0 1.0; 
        0.0 0.0 0.0 -1.0]) < 1.0e-8
end