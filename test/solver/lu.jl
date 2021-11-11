@testset "LU solver" begin
    n = 20
    m = 10
    A = rand(n, n)
    X = rand(n, m)
    B = rand(n, m)
    x = rand(n)
    b = rand(n)

    solver = RobotDojo.empty_solver(A) 
    @test solver isa RobotDojo.LinearSolver

    solver = RobotDojo.lu_solver(A)
    RobotDojo.linear_solve!(solver, X, A, B)
    @test norm(A * X - B, Inf) < 1.0e-10

    solver = RobotDojo.lu_solver(A)
    RobotDojo.linear_solve!(solver, x, A, b)
    @test norm(A * x - b, Inf) < 1.0e-10


end