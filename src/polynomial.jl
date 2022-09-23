#############################################################################
#############################################################################
#
# This file defines the abstract polynomial type along with several generic operations 
#                                                                               
#############################################################################
#############################################################################
abstract type Polynomial end

"""
Construct a polynomial of the form x^p-x.
"""
cyclotonic_polynomial(polyType::Type{<:Polynomial},p::Int) = polyType([Term(1,p), Term(-1,0)])

"""
Construct a polynomial of the form x-n.
"""
linear_monic_polynomial(polyType::Type{<:Polynomial},n::Integer) = polyType([Term(1,1), Term(-n,0)])

"""
Construct a polynomial of the form x.
"""
x_poly(polyType::Type{<:Polynomial}) = polyType(Term(1,1))

"""
Creates the zero polynomial.
"""
zero(polyType::Type{<:Polynomial}) = polyType()
zero(p::Polynomial) = zero(typeof(p))

"""
Creates the unit polynomial.
"""
one(polyType::Type{<:Polynomial}) = polyType(one(Term))
one(p::Polynomial) = one(typeof(p))

"""
Generates a random polynomial.
"""
function rand(polyType::Type{<:Polynomial} ; 
                degree::Int = -1, 
                terms::Int = -1, 
                max_coeff::Int = 100, 
                mean_degree::Float64 = 5.0,
                prob_term::Float64  = 0.7,
                monic = false,
                condition = (p)->true)
        
    while true 
        _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
        _terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
        degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
        coeffs = rand(1:max_coeff,_terms+1)
        monic && (coeffs[end] = 1)
        p = polyType( [Term(coeffs[i],degrees[i]) for i in 1:length(degrees)] )
        condition(p) && return p
    end
end



##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the non-zero terms of the polynomial. This implements the iteration interface.
"""
iterate(p::Polynomial, state=1) = iterate(p.terms, state)

##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::Polynomial) = length(p.terms) 

"""
Returns the coefficients of the polynomial.
"""
coeffs(p::Polynomial)::Vector{Integer} = [t.coeff for t in p]

"""
The degree of the polynomial.
"""
degree(p::Polynomial)::Integer = leading(p).degree 

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::Polynomial)::Integer = euclid_alg(coeffs(p))

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(f::Polynomial, x::T) where T <: Number = sum(evaluate(t,x) for t in f)

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

"""
Copy a polynomial
"""

copyTerms(p::Polynomial) = typeof(p)(copy(p.terms),sorting=false,filtering=false)


"""
The negative of a polynomial.
"""
-(p::Polynomial) = typeof(p)(map((pt)->-pt, p.terms))

"""
Create a new polynomial which is the derivative of the polynomial.
"""
function derivative(p::Polynomial)::Polynomial
    #make sure the output type matches the input type
    der_p = zero(p)
    for term in p
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p,der_term)
    end
    return trim!(der_p)
end


"""
The prim part (multiply a polynomial by the inverse of its content).
"""
# function prim_part(p::Polynomial)
#     c = content(p)
#     d = last(p.terms).degree 
#     reduced = copyTerms(p)
#     map!(t->typeof(t)(t.coeff, t.degree-d),reduced.terms,reduced.terms)
#     reduced = reduced ÷ c
#     return reduced
# end

prim_part(p::Polynomial) = p ÷ content(p)


"""
A square free polynomial.
"""
square_free(p::Polynomial, prime::Integer)::Polynomial = (p ÷ gcd(p,derivative(p),prime))(prime)

#################################
# Queries about two polynomials #
#################################

"""
Check if two polynomials are the same
"""
==(p1::Polynomial, p2::Polynomial)::Bool = p1.terms == p2.terms


"""
Check if a polynomial is equal to 0.
"""
#Note that in principle there is a problem here. E.g The polynomial 3 will return true to equalling the integer 2.
==(p::Polynomial, n::T) where T <: Real = iszero(p) == iszero(n)

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Subtraction of two polynomials.
"""
-(p1::Polynomial, p2::Polynomial)::Polynomial = p1 + (-p2)


"""
Multiplication of polynomial and term.
"""
function mult!(p::Polynomial,t::Term)::Polynomial
    if iszero(t)
        empty!(p.terms)
    else
        for i in 1:length(p)
            p.terms[i] *= t
        end
    end
    return p
end

*(p::Polynomial,t::Term) = mult!(copyTerms(p),t)
*(t::Term,p::Polynomial) = p*t

"""
Multiplication of polynomial and an integer.
"""
*(n::Integer, p::Polynomial)::Polynomial = p*Term(n,0)
*(p::Polynomial, n::Integer)::Polynomial = n*p

"""
Integer division of a polynomial by an integer.

Warning this may not make sense if n does not divide all the coefficients of p.
"""
÷(p::Polynomial, n::Integer) = (prime)->typeof(p)(map((pt)->((pt ÷ n)(prime)), p.terms))

"""
Take the mod of a polynomial with an integer.
"""
function mod!(f::Polynomial, p::Integer)::Polynomial
    for i in 1:length(f.terms)
        f.terms[i] = mod(f.terms[i], p)
    end
    return trim!(f)
end

mod(f::Polynomial, p::Integer) = mod!(copyTerms(f),p)

function smod(f::Polynomial, p::Integer)::Polynomial
    f_out = copyTerms(f)
    for i in 1:length(f_out.terms)
        f_out.terms[i] = smod(f_out.terms[i], p)
    end
    return trim!(f_out)
end

"""
Power of a polynomial mod prime.
"""
function pow_mod(p::Polynomial, n::Integer, prime::Integer)
    out = PolynomialModP(p,prime)
    return (out^n).polynomial
end