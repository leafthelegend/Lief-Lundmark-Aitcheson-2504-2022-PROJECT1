#############################################################################
#############################################################################
#
# This file contains units tests for polynomial factorization
#                                                                               
#############################################################################
#############################################################################

using Primes

function factor_test_poly(T::Type{<:Polynomial};N::Int = 10, seed::Int = 0, primes::Vector{Int} = [5,17,19])
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            p = rand(T)
            factorization = factor(p, prime)
            pr = mod(expand_factorization(factorization),prime)
            @assert mod(p-pr,prime) == 0 
        end
    end

    println("\nfactor_test_poly - PASSED")
end


"""
Test factorization of polynomials.
"""
function factor_test_poly1(T::Type{<:Polynomial};N::Int = 10, seed::Int = 0, N_primes::Int=100,lowest_prime::Int = 200,mean_degree::Float64 = 10.,max_coeff::Integer = 100)
    
    if length(T.parameters) > 0
        max_coeff = T.parameters[1](max_coeff)
    end
    
    Random.seed!(seed)
    prime = lowest_prime
    for _ in 1:N_primes
        prime = Primes.nextprime(prime+1)
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            p = rand(T,mean_degree=mean_degree,max_coeff = max_coeff)
            factorization = factor(p, prime)
            pr = mod(expand_factorization(factorization),prime)
            @assert mod(p-pr,prime) == 0 
            print(".")
        end
    end

    println("\nfactor_test_poly - PASSED")
end

function get_reducible_poly(T::Type{<:Polynomial},degree::Int,prime::Int,max_coeff::Integer)::PolynomialModP
    factors = []
    total_degree = 0
    n = 1
    while total_degree<degree
        push!(factors,rand_modp(T,prime,degree=n,max_coeff = max_coeff))
        total_degree += n
        n += 1
    end
    p=factors[1]
    for i in 2:length(factors)
        p *= factors[i]
    end
    return p
end

function factor_test_poly2(T::Type{<:Polynomial};N::Int = 10, seed::Int = 0, N_primes::Int=100,lowest_prime::Int = 19,degree::Int = 10,max_coeff::Integer = 10000)
    
    if length(T.parameters) > 0
        max_coeff = T.parameters[1](max_coeff)
    end
    
    Random.seed!(seed)
    prime = lowest_prime
    for _ in 1:N_primes
        prime = Primes.nextprime(prime+1)
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            p = get_reducible_poly(T,degree,prime,max_coeff)
            println(p)
            factorization = factor(p)
            println(factorization)
            pr = expand_factorization(factorization)
            @assert p == pr
            print(".")
        end
    end

    println("\nfactor_test_poly2 - PASSED")
end