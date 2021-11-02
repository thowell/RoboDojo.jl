@testset "Solver: Random QP (LU)" begin
    """
        minimize   x' P x + q' x
        subject to    x >= 0
    """

    n = 100

    _P = rand(n)
    P = Diagonal(_P)
    q = rand(n)
    θ = [_P; q]
    z = ones(2 * n)
    r = zeros(2 * n)
    rz = zeros(2 * n, 2 * n)
    rθ = zeros(2 * n, 2 * n)
    κ = [1.0]

    # residual
    function _r!(r, z, θ, κ)
        x = z[1:100]
        y = z[101:200]
        P = θ[1:100]
        q = θ[101:200]

        r[1:100] = 2.0 * Diagonal(P) * x + q - y
        r[101:200] = Diagonal(x) * y .- κ[1]
        nothing
    end

    @variables r_sym[1:200]
    @variables z_sym[1:200]
    @variables θ_sym[1:200]
    @variables κ_sym[1:1]

    parallel = Symbolics.SerialForm()
    _r!(r_sym, z_sym, θ_sym, κ_sym)
    rf! = eval(Symbolics.build_function(r_sym, z_sym, θ_sym, κ_sym,
        parallel = parallel)[2])
    rz_exp = Symbolics.jacobian(r_sym, z_sym)
    rθ_exp = Symbolics.jacobian(r_sym, θ_sym)
    rz_sp = similar(rz_exp, Float64)
    rθ_sp = similar(rθ_exp, Float64)
    rzf! = eval(Symbolics.build_function(rz_exp, z_sym, θ_sym,
        parallel = parallel)[2])
    rθf! = eval(Symbolics.build_function(rθ_exp, z_sym, θ_sym,
        parallel = parallel)[2])

    # options
    opts = RobotDojo.InteriorPointOptions(diff_sol = true, verbose=false)
    idx = RobotDojo.IndicesOptimization(
            2n, 2n, 
            [collect(1:n), collect(n .+ (1:n))],
            [collect(1:n), collect(n .+ (1:n))], 
            Vector{Vector{Int}}[], 
            Vector{Vector{Int}}[], 
            collect(1:n), 
            collect(n .+ (1:n)), 
            collect(1:0),
            Vector{Int}[],
            collect(n .+ (1:n)))

    # solver
    ip = RobotDojo.interior_point(z, θ,
        idx=idx,
        r! = rf!, rz! = rzf!, rθ! = rθf!,
        rz = rz_sp,
        rθ = rθ_sp,
        opts = opts)

    # solve
    status = interior_point_solve!(ip)

    # test
    @test status
    @test norm(ip.r, Inf) < opts.r_tol
    @test all(ip.z[idx.ortz[1]] .> 0.0)
    @test all(ip.z[idx.ortz[2]] .> 0.0)
    @test ip.κ[1] < opts.κ_tol
    @test norm(ip.δz, 1) != 0.0
end
