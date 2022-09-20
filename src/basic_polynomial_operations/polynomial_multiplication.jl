#############################################################################
#############################################################################
#
# This file implements polynomial multiplication 
#                                                                               
#############################################################################
#############################################################################

"""
Multiply two polynomials.
"""
function *(p1::Polynomial, p2::Polynomial)::Polynomial
    p_out = typeof(p1)()
    for t in p1
        new_summand = (t * p2)
        p_out = add!(p_out, new_summand)
    end
    return p_out
end

"""
Divide and conquer
"""
function split(p::Polynomial)::Tuple{Polynomial,Polynomial}
    n = length(p)
    if n == 1
        return (p[1],typeof(p)())
    end
    n2 = n÷2
    p1 = typeof(p)(p.terms[1:n2])
    p2 = typeof(p)(p.terms[n2+1:end])
    return (p1,p2)
end

function fast_multiply(p1::Polynomial, p2::Polynomial)::Polynomial
    if iszero(p1) || iszero(p2)
        return typeof(p1)()
    end
    if length(p1) == 1
        return p1.terms[1] * p2
    end
    if length(p2) == 1
        return p2.terms[1] * p1
    end
    p1_1, p1_2 = split(p1)
    p2_1, p2_2 = split(p2)
    return add!(add!(fast_multiply(p1_1,p2_1), fast_multiply(p1_1,p2_2)), add!(fast_multiply(p1_2,p2_1), fast_multiply(p1_2,p2_2)))
end

# *(p1::Polynomial, p2::Polynomial) = fast_multiply(p1,p2)

"""
Power of a polynomial.
"""
function ^(p::Polynomial, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end

