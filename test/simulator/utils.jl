@testset "Simulator: utilities" begin
    # dimensions
    nz = RoboDojo.num_var(biped) 
    nθ = RoboDojo.num_data(biped, nf=length(friction_coefficients(biped)))

    @test nz == 57 
    @test nθ == 33 

    # trajectory 
    traj = RoboDojo.Trajectory(biped, 5) 
    grad = RoboDojo.GradientTrajectory(biped, 5)
    @test length(traj.q) == 7
    @test length(traj.v) == 6
    @test length(traj.u) == 5
    traj.q[1] .= rand(biped.nq)
    RoboDojo.reset!(traj)
    @test all(traj.q[1] .== 0.0)
    @test length(grad.∂q3∂q1) == 5
    grad.∂q3∂q1[1] .= rand(biped.nq, biped.nq)
    RoboDojo.reset!(grad)
    @test all(grad.∂q3∂q1[1] .== 0.0)

    # disturbances 
    dist = RoboDojo.empty_disturbances(biped) 
    @test length(dist.w) == biped.nw
    @test length(RoboDojo.disturbances(dist, traj, 1)) == biped.nw

    # policy 
    p = RoboDojo.empty_policy(biped)
    @test length(RoboDojo.policy(p, traj, 1)) == biped.nu

    # indices
    
    idx_z = RoboDojo.indices_z(biped)
    @test idx_z.q == collect(1:9)
    @test idx_z.γ == collect(10:17) 
    @test idx_z.sγ == collect(18:25)
    @test idx_z.ψ == collect(26:33)
    @test idx_z.b == collect(34:41)
    @test idx_z.sψ == collect(42:49)
    @test idx_z.sb == collect(50:57)

    idx_θ = RoboDojo.indices_θ(biped, nf=length(friction_coefficients(biped))) 
    @test idx_θ.q1 == collect(1:9)
    @test idx_θ.q2 == collect(10:18)
    @test idx_θ.u == collect(19:24)
    @test idx_θ.w == collect(25:24)
    @test idx_θ.f == collect(25:32)
    @test idx_θ.h == collect(33:33)

    idx = RoboDojo.IndicesOptimization() 
    @test idx.nz == 0 
    @test idx.nΔ == 0 

    idx_opt = RoboDojo.indices_optimization(biped)
    @test idx_opt.nz == nz
    @test idx_opt.nΔ == nz
    @test idx_opt.ortz == [[10, 11, 12, 13, 14, 15, 16, 17], [18, 19, 20, 21, 22, 23, 24, 25]]
    @test idx_opt.ortΔ == [[10, 11, 12, 13, 14, 15, 16, 17], [18, 19, 20, 21, 22, 23, 24, 25]]
    @test idx_opt.socz == [[[26, 34], [42, 50]],
                            [[27, 35], [43, 51]],
                            [[28, 36], [44, 52]],
                            [[29, 37], [45, 53]],
                            [[30, 38], [46, 54]],
                            [[31, 39], [47, 55]],
                            [[32, 40], [48, 56]],
                            [[33, 41], [49, 57]]]
    @test idx_opt.socΔ == [[[26, 34], [42, 50]],
                            [[27, 35], [43, 51]],
                            [[28, 36], [44, 52]],
                            [[29, 37], [45, 53]],
                            [[30, 38], [46, 54]],
                            [[31, 39], [47, 55]],
                            [[32, 40], [48, 56]],
                            [[33, 41], [49, 57]]]
    @test idx_opt.equr == collect(1:33)
    @test idx_opt.ortr == collect(34:41)
    @test idx_opt.socr == collect(42:57)
    @test idx_opt.socri == [[42, 43],
                            [44, 45],
                            [46, 47],
                            [48, 49],
                            [50, 51],
                            [52, 53],
                            [54, 55],
                            [56, 57]]                
    @test idx_opt.bil == collect(34:57)

    # initialize 
    z0 = zeros(nz) 
    θ0 = zeros(nθ)
    q0 = rand(biped.nq)
    u0 = rand(biped.nu) 
    w0 = rand(biped.nw) 
    f0 = friction_coefficients(biped) 
    h0 = [0.235] 
    RoboDojo.initialize_z!(z0, biped, idx_z, q0)
    @test norm(z0 - [q0; 1.0 * ones(3 * biped.nc); 0.1 * ones(biped.nc); 1.0 * ones(biped.nc); 0.1 * ones(biped.nc)]) < 1.0e-8

    RoboDojo.initialize_θ!(θ0, biped, idx_θ, q0, q0, u0, w0, f0, h0)
    @test norm(θ0 - [q0; q0; u0; w0; f0; h0]) < 1.0e-8

    # floating base dimensions
    @test RoboDojo.floating_base_dim(biped) == 3
end