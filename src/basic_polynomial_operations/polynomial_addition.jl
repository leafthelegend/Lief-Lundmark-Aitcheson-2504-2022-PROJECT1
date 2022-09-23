#############################################################################
#############################################################################
#
# This file implements polynomial addition 
#                                                                               
#############################################################################
#############################################################################

"""
Add a dense polynomial and a term (in place).
"""
function add!(p::PolynomialDense, t::Term)
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
Add a dense polynomial and a term (not in place).
"""
+(p::PolynomialDense, t::Term) = add!(copyTerms(p), t)



"""
Add a sparse polynomial and a term (in place).
"""
function add!(p::PolynomialSparse, t::Term)
    if isempty(p.terms)
        push!(p.terms,t)
        return p
    end
    i = findfirst(term -> term.degree <= t.degree, p.terms)
    if isnothing(i)
        push!(p.terms,t)
    elseif t.degree == p.terms[i].degree
        p.terms[i] += t
        if iszero(p.terms[i])
            deleteat!(p.terms,i)
        end
    else
        insert!(p.terms,i,t)
    end
    return p
end

"""
Add a sparse polynomial and a term (not in place).
"""
+(p::PolynomialSparse, t::Term) = add!(copyTerms(p), t)

+(t::Term, p::Polynomial) = p + t

"""
Add two sparse polynomials (in place).
"""
function add!(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    t1, t2 = length(p1),length(p2) #initialise pointers to the terms of p1 and p2
    while t2 > 0
        if t1 == 0
            prepend!(p1.terms, p2.terms[1:t2])
            break
        end
        if p1.terms[t1].degree > p2.terms[t2].degree
            insert!(p1.terms,t1+1,p2.terms[t2])
            t2 -= 1
        elseif p1.terms[t1].degree < p2.terms[t2].degree
            t1 -= 1
        else
            p1.terms[t1] += p2.terms[t2]
            if iszero(p1.terms[t1])
                deleteat!(p1.terms,t1)
            end
            t1 -= 1
            t2 -= 1
        end
    end
    return p1
end

"""
Add two dense polynomials (in place)
"""
function add!(p1::PolynomialDense,p2::PolynomialDense)::PolynomialDense
    for t in p2
        p1+=t
    end
    return p1
end

"""
Add two polynomials (not in place).
"""
+(p1::Polynomial, p2::Polynomial) = add!(copyTerms(p1), p2)

"""
Add a polynomial and an integer.
"""
+(p::PolynomialDense, n::Int) = p + Term(n,0)
+(n::Int, p::PolynomialDense) = p + Term(n,0)
-(p::PolynomialDense, n::Int) = p + Term(-n,0)
-(n::Int, p::PolynomialDense) = -p + Term(n,0)

#ensure integer is converted to correct type for addition with PolynomialSparse

+(p::PolynomialSparse{T}, n::Integer) where {T<:Integer} = p + Term{T}(n,0)
+(n::Integer, p::PolynomialSparse{T}) where {T<:Integer} = p + Term{T}(n,0)
-(p::PolynomialSparse{T}, n::Integer) where {T<:Integer} = p + Term{T}(-n,0)
-(n::Integer, p::PolynomialSparse{T}) where {T<:Integer} = -p + Term{T}(n,0)