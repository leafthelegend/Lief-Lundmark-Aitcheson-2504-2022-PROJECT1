#############################################################################
#############################################################################
#
# This file contains unit tests for PolynomialModP operations
#                                                                               
#############################################################################
#############################################################################
using Primes

"""
Test product of polynomials.
"""
function prod_test_polyModP(T::Type{<:Polynomial};N::Int = 10^3, N_prods::Int = 20, seed::Int = 0, maxPrime::Int = 50, verbose::Bool = true)
    Random.seed!(seed)
    primes = Primes.primes(maxPrime)

    for _ in 1:N
        prime = Random.rand(primes)
        p1 = PolynomialModP(rand(T),prime)
        p2 = PolynomialModP(rand(T),prime)
        prod = p1*p2
        @assert leading(prod.polynomial) == mod((leading(p1.polynomial)*leading(p2.polynomial)),prime)
    end

    for _ in 1:N
        prime = Random.rand(primes)
        if length(T.parameters)>0
            p_base = T(Term{T.parameters[1]}(1,0))
        else
            p_base = T(Term(1,0))
        end
        p_base = PolynomialModP(p_base,prime)
        for _ in 1:N_prods
            p = rand(T)
            prod = p_base*PolynomialModP(p,prime)
            #we must also handle the case where the leading term is zero mod p, so the first part of this assertion will fail
            @assert leading(prod.polynomial) == mod((leading(p_base.polynomial)*leading(p)),prime) || iszero(mod((leading(p_base.polynomial)*leading(p)),prime))
            p_base = prod
        end
    end
    verbose && println("prod_test_polyModP - PASSED")
    nothing
end

"""
Test division of polynomials modulo p.
"""
function division_test_polyModP(T::Type{<:Polynomial}; N::Int = 10^4, seed::Int = 0, maxPrime::Int=50, verbose::Bool=true)
    Random.seed!(seed)
    primes = Primes.primes(maxPrime)
    
    for _ in 1:N
        prime = Random.rand(primes)
        p1 = PolynomialModP(rand(T),prime)
        p2 = PolynomialModP(rand(T),prime)
        p_prod = p1*p2
        q, r = PolynomialModP(T(),prime), PolynomialModP(T(),prime)
        try
            q, r = divide(p_prod, p2)
            if (q.polynomial, r.polynomial) == (nothing,nothing)
                verbose && println("Unlucky prime: numerator is $p1, denominator is $p2")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert p2 == 0
            else
                throw(e)
            end
        end
        @assert iszero(q*p2+r - p_prod)
    end
    verbose && println("division_test_polyModP - PASSED")
end

"""
Test factorization of polynomials mod p.
"""
function factor_test_polyModP(T::Type{<:Polynomial};N::Int = 10, seed::Int = 0, primes::Vector{Int} = [5,17,19], verbose::Bool=true)
    Random.seed!(seed)
    for prime in primes
        verbose && print("\ndoing prime = $prime \t")
        for _ in 1:N
            verbose && print(".")
            p = PolynomialModP(rand(T),prime)
            factorization = factor(p)
            pr = expand_factorization(factorization)
            @assert p == pr
        end
    end

    println("\nfactor_test_polyModP - PASSED")
end