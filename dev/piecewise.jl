using Symbolics
using BenchmarkTools

@variables x[1:1]

m = max.(x, 0.0)
function height(x) 
    [Symbolics.IfElse.ifelse(x[1] >= 0.0, 
        Symbolics.IfElse.ifelse(x[1] >= 2.0, 2.0 * x[1], 1.0 * x[1]), 
        0.0)]
end

m = height(x)

dm = Symbolics.jacobian(m, x)

m_func = eval(Symbolics.build_function(m, x)[1])
m_func! = eval(Symbolics.build_function(m, x)[2])
dm_func = eval(Symbolics.build_function(dm, x)[1])
dm_func! = eval(Symbolics.build_function(dm, x)[2])

a = [-1.0]

m_cache = zeros(length(m))
dm_cache = zeros(length(dm))

@benchmark m_func($a)
@benchmark m_func!($m_cache, $a)
@benchmark dm_func($a)
@benchmark dm_func!($dm_cache, $a)


