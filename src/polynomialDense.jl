#############################################################################
#############################################################################
#
# This file defines the dense polynomial type with several operations 
#                                                                               
#############################################################################
#############################################################################

####################################
# PolynomialDense type and construction #
####################################

"""
A PolynomialDense type - designed to be for polynomials with integer coefficients.
"""
struct PolynomialDense <: Polynomial

    #A zero packed vector of terms
    #Terms are assumed to be in order with first term having degree 0, second degree 1, and so fourth
    #until the degree of the polynomial. The leading term (i.e. last) is assumed to be non-zero except 
    #for the zero polynomial where the vector is of length 1.
    #Note: at positions where the coefficient is 0, the power of the term is also 0 (this is how the Term type is designed)
    terms::Vector{Term}   
    
    #Inner constructor of 0 polynomial
    PolynomialDense() = new([zero(Term)])

    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialDense(vt::Vector{Term})

        #Filter the vector so that there is not more than a single zero term
        vt = filter((t)->!iszero(t), vt)
        if isempty(vt)
            vt = [zero(Term)]
        end

        max_degree = maximum((t)->t.degree, vt)
        terms = [zero(Term) for i in 0:max_degree] #First set all terms with zeros

        #now update based on the input terms
        for t in vt
            terms[t.degree + 1] = t #+1 accounts for 1-indexing
        end
        return new(terms)
    end
end

"""
This function maintains the invariant of the PolynomialDense type so that there are no zero terms beyond the highest
non-zero term.
"""
function trim!(p::PolynomialDense)::PolynomialDense
    i = length(p.terms)
    while i > 1
        if iszero(p.terms[i])
            pop!(p.terms)
        else
            break
        end
        i -= 1
    end
    return p
end

"""
Check if the dense polynomial is zero.
"""
iszero(p::PolynomialDense)::Bool = p.terms == [Term(0,0)]

"""
Construct a dense polynomial with a single term.
"""
PolynomialDense(t::Term) = PolynomialDense([t])

###########
# Display #
###########

"""
Show a dense polynomial.
"""
function show(io::IO, p::PolynomialDense)
    if iszero(p)
        print(io,"0")
    else
        n = length(p.terms)
        for (i,t) in Iterators.reverse(enumerate(p.terms))
            if !iszero(t)
                print(io, i == n ? t : "$(t.coeff>0 ? "+" : "")$(t)")
            end
        end
    end
end

##############################
# Queries about a dense polynomial #
##############################

"""
The leading term of the polynomial.
"""
leading(p::PolynomialDense)::Term = isempty(p.terms) ? zero(Term) : last(p.terms) 

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the polynomial.
"""
#Note that ideally this would throw and error if pushing another term of degree that is already in the polynomial
function push!(p::PolynomialDense, t::Term) 
    if t.degree <= degree(p)
        p.terms[t.degree + 1] = t
    else
        append!(p.terms, zeros(Term, t.degree - degree(p)-1))
        push!(p.terms, t)
    end
    return p        
end

"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.
"""
function pop!(p::PolynomialDense)::Term 
    popped_term = pop!(p.terms) #last element popped is leading coefficient

    while !isempty(p.terms) && iszero(last(p.terms))
        pop!(p.terms)
    end

    if isempty(p.terms)
        push!(p.terms, zero(Term))
    end

    return popped_term
end

