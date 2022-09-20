#############################################################################
#############################################################################
#
# This file implements polynomial GCD 
#                                                                               
#############################################################################
#############################################################################

"""
The extended euclid algorithm for polynomials modulo prime. s and t will take the type of a. g may take the type of either a or b.
"""
function extended_euclid_alg(a::Polynomial, b::Polynomial, prime::Integer)
    old_r, r = mod(a, prime), mod(b, prime)
    old_s, s = one(a), zero(a)
    old_t, t = zero(a), one(a)

    while !iszero(mod(r,prime))
        q = first(divide(old_r, r)(prime))
        old_r, r = r, mod(old_r - q*r, prime)
        old_s, s = s, mod(old_s - q*s, prime)
        old_t, t = t, mod(old_t - q*t, prime)
    end
    g, s, t = old_r, old_s, old_t

    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t  
end

"""
The GCD of two polynomials modulo prime.
"""
gcd(a::Polynomial, b::Polynomial, prime::Integer) = extended_euclid_alg(a,b,prime) |> first