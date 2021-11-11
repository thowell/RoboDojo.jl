@testset "Utilities" begin
    @test norm(RobotDojo.rotation_matrix(0.0) - I(2)) < 1.0e-8
    @test norm(RobotDojo.vector_rotation_matrix([0.0; 1.0]) - I(2)) < 1.0e-8
    @test RobotDojo.scn(1.0, digits=2) == " 1.00e+0"
    @test RobotDojo.scn(12.345, digits=2) == " 1.23e+1"
end