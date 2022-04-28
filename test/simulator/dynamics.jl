@testset "Simulator: dynamics" begin
    dynamics_model = Simulator(halfcheetah, 1, h=h)
    dynamics_model.ip.opts.r_tol = 1e-7
    dynamics_model.ip.opts.κ_tol = 1e-5
    dynamics_model.ip.opts.undercut = 10.0
    nq = dynamics_model.model.nq
    nx = 2nq
    nu = dynamics_model.model.nu
    nw = dynamics_model.model.nw

    y = zeros(nx)
    dx = zeros(nx,nx)
    du = zeros(nx,nu)
    x = RoboDojo.nominal_state(dynamics_model.model) + 0.1*ones(nx)
    u = ones(nu)
    w = ones(nw)

    status = RoboDojo.dynamics(dynamics_model, y, x, u, w)
    @test status
    status = RoboDojo.dynamics_jacobian_state(dynamics_model, dx, x, u, w)
    @test status
    status = RoboDojo.dynamics_jacobian_input(dynamics_model, du, x, u, w)
    @test status

    function out_of_place_explicit_dynamics(dynamics_model, x, u, w)
        nx = 2dynamics_model.model.nq
        y = zeros(nx)
        RoboDojo.dynamics(dynamics_model, y, x, u, w)
        return y
    end

    dx0 = FiniteDiff.finite_difference_jacobian(x -> out_of_place_explicit_dynamics(dynamics_model, x, u, w), x)
    du0 = FiniteDiff.finite_difference_jacobian(u -> out_of_place_explicit_dynamics(dynamics_model, x, u, w), u)
    @test norm(dx - dx0, Inf) < 1e-4
    @test norm(du - du0, Inf) < 1e-4

end
