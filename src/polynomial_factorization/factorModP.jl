#############################################################################
#############################################################################
#
# This file implements factorization 
#                                                                               
#############################################################################
#############################################################################

"""
Factors a polynomial over the field Z_p.

Returns a vector of tuples of (irreducible polynomials (mod p), multiplicity) such that their product of the list (mod p) is f. Irreducibles are fixed points on the function factor.
"""
function factor(f::PolynomialModP)::Vector{Tuple{PolynomialModP,Integer}}
    #Cantor Zassenhaus factorization
    degree(f) ≤ 1 && return [(f,1)]

    # make f primitive
    ff = prim_part(f)
    # @show "after prim:", ff

     # make f square-free
    squares_poly = gcd(f, derivative(ff))
    # @show "squares_poly:", squares_poly
    ff = ff ÷ squares_poly
    # @show "after square free:", ff

    # make f monic
    old_coeff = leading(ff).coeff
    ff = ff ÷ old_coeff
    # @show "after monic:", ff

    dds = dd_factor(ff)

    ret_val = Tuple{PolynomialModP,Integer}[]
    for (k,dd) in enumerate(dds)
        sp = dd_split(dd, k)
        sp = map((p)->p ÷ leading(p).coeff,sp) #makes the polynomials inside the list sp, monic
        for mp in sp
            push!(ret_val, (mp, multiplicity(f,mp)) )
        end
    end

    # #account for factors of x^n which were removed by the square-free reduction
    # if multiplicity(f,x_poly(f)) > 0
    #     push!(ret_val, (x_poly(f), multiplicity(f,x_poly(f))) )
    # end

    #Append the leading coefficient as well
    push!(ret_val, (leading(f).coeff* one(f), 1) )

    return ret_val
end

"""
Compute the number of times g divides f
"""
function multiplicity(f::PolynomialModP, g::PolynomialModP)::Integer
    degree(gcd(f, g)) == 0 && return 0
    return 1 + multiplicity(f ÷ g, g)
end


"""
Distinct degree factorization.

Given a square free polynomial `f` returns a list, `g` such that `g[k]` is a product of irreducible polynomials of degree `k` for `k` in 1,...,degree(f) ÷ 2, such that the product of the list (mod `prime`) is equal to `f` (mod `prime`).
"""
function dd_factor(f::PolynomialModP)::Array{PolynomialModP}
    x = x_poly(f)
    w = copyTerms(x)
    g = Array{typeof(f)}(undef,degree(f)) #Array of polynomials indexed by degree

    #Looping over degrees
    for k in 1:degree(f)
        w = rem(w^f.p, f)
        g[k] = gcd(w - x, f) 
        f = f ÷ g[k]
    end

    #edge case for final factor
    f != one(f) && push!(g,f)
    
    return g
end

"""
Distinct degree split.

Returns a list of irreducible polynomials of degree `d` so that the product of that list (mod prime) is the polynomial `f`.
"""
function dd_split(f::PolynomialModP, d::Integer)::Vector{PolynomialModP}
    degree(f) == d && return [f]
    degree(f) == 0 && return []
    w = rand_modp(typeof(f.polynomial), f.p, degree = d, monic = true)
    n_power = (f.p^d-1) ÷ 2
    wp = w^n_power - one(f)
    g = gcd(wp, f)
    ḡ = f ÷ g
    return vcat(dd_split(g, d), dd_split(ḡ, d) )
end

"""
Expand a factorization.
"""
function expand_factorization(factorization::Vector{<:Tuple{Union{Polynomial,PolynomialModP},Integer}})::Union{Polynomial,PolynomialModP}
    length(factorization) == 1 && return first(factorization[1])^last(factorization[1])
    return *([first(tt)^last(tt) for tt in factorization]...)
end

function factor(f::Polynomial,prime::Int)::Vector{Tuple{Polynomial,Integer}}
    factorisation = factor(PolynomialModP(f, prime))
    return [(p.polynomial, m) for (p,m) in factorisation]
end