function cone_product(z, s)
    n = length(z)
    [z' * s; z[1] * view(s, 2:n) + s[1] * view(z, 2:n)]
end

function cone_product!(a, z, s; reset=false) 
    reset && fill!(a, 0.0)
    n = length(a)
    a[1] += z' * s
    for i = 2:n
        a[i] += z[1] * s[i]
        a[i] += s[1] * z[i]
    end
    return nothing
end