"""
Power of a polynomial.
"""
function ^(p::Union{Polynomial,PolynomialModP}, n::Int)
    max_power = floor(Int, log2(n))
    powers::Vector{typeof(p)} = [p] #p^1, p^2, p^4, p^8, p^16, ...
    for i in 1:max_power
        powers = push!(powers, powers[end]*powers[end])
    end
    out = one(p)
    for i in 0:max_power
        if n & (1 << i) != 0 #if n has a 1 in the ith bit
            out*=powers[i+1]
        end
    end
    return out
end

function pow_mod(p::Polynomial, n::Integer, prime::Integer)
    out = PolynomialModP(p,prime)
    return (out^n).polynomial
end