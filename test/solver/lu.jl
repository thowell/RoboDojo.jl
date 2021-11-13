@testset "LU solver" begin
    n = 20
    m = 10
    A = rand(n, n)
    X = rand(n, m)
    B = rand(n, m)
    x = rand(n)
    b = rand(n)

    solver = RoboDojo.empty_solver(A) 
    @test solver isa RoboDojo.LinearSolver

    solver = RoboDojo.lu_solver(A)
    RoboDojo.linear_solve!(solver, X, A, B)
    @test norm(A * X - B, Inf) < 1.0e-10

    solver = RoboDojo.lu_solver(A)
    RoboDojo.linear_solve!(solver, x, A, b)
    @test norm(A * x - b, Inf) < 1.0e-10


end