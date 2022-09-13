#############################################################################
#############################################################################
#
# This file implements polynomial addition 
#                                                                               
#############################################################################
#############################################################################

"""
Add a dense polynomial and a term.
"""
function +(p::PolynomialDense, t::Term)
    p = deepcopy(p)
    if t.degree > degree(p)
        push!(p, t)
    else
        if !iszero(p.terms[t.degree + 1]) #+1 is due to indexing
            p.terms[t.degree + 1] += t
        else
            p.terms[t.degree + 1] = t
        end
    end

    return trim!(p)
end

"""
Add a sparse polynomial and a term.
"""
function +(p::PolynomialSparse, t::Term)
    p = deepcopy(p)
    if isempty(p.terms)
        push!(p.terms,t)
        return p
    end
    i = findfirst(term -> term.degree <= t.degree, p.terms)
    if isnothing(i)
        push!(p.terms,t)
    elseif t.degree == p.terms[i].degree
        p.terms[i] += t
    else
        insert!(p.terms,i,t)
    end
    return trim!(p)
end

+(t::Term, p::Polynomial) = p + t

"""
Add two polynomials.
"""
function +(p1::Polynomial, p2::Polynomial)::Polynomial
    p = deepcopy(p1)
    for t in p2
        p += t
    end
    return p
end

"""
Add a polynomial and an integer.
"""
+(p::Polynomial, n::Int) = p + Term(n,0)
+(n::Int, p::Polynomial) = p + Term(n,0)
-(p::Polynomial, n::Int) = p + Term(-n,0)
-(n::Int, p::Polynomial) = -p + Term(n,0)

+(p::PolynomialSparse{T}, n::T) where {T<:Integer} = p + Term{T}(n,0)
+(n::T, p::PolynomialSparse{T}) where {T<:Integer} = p + Term{T}(n,0)
-(p::PolynomialSparse{T}, n::T) where {T<:Integer} = p + Term{T}(-n,0)
-(n::T, p::PolynomialSparse{T}) where {T<:Integer} = -p + Term{T}(n,0)