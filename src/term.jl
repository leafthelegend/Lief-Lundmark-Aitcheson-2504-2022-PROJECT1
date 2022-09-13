#############################################################################
#############################################################################
#
# This file defines the Term type with several operations 
#                                                                               
#############################################################################
#############################################################################

##############################
# Term type and construction #
##############################

"""
A term.
"""
struct Term{T<:Integer}  #structs are immutable by default
    coeff::T
    degree::Int
    function Term{T}(coeff::Integer, degree::Int) where {T<:Integer}
        degree < 0 && error("Degree must be non-negative")
        coeff != 0 ? new{T}(coeff,degree) : new{T}(coeff,0)
    end
end

Term(coeff::T,degree::Int) where {T<:Integer} = Term{T}(coeff,degree)

"""
Creates the zero term.
"""
zero(::Type{<:Term{T}}) where {T} = Term{T}(0,0)
zero(::Type{Term}) = Term(0,0) #default to Int


"""
Creates the unit term.
"""
one(::Type{<:Term{T}}) where {T} = Term{T}(1,0)
one(::Type{Term}) = Term(1,0) #default to Int

###########
# Display #
###########

"""
Show a term.
"""
function show(io::IO, t::Term, cdot::Bool = false)
    iszero(t.coeff) && return(print(io,"0"))
    iszero(t.degree) && return(print(io,t.coeff))
    print(io,t.coeff == 1 ? "" : "$(t.coeff)$(cdot ? "⋅" : "")",(t.degree > 1 ? "x^$(t.degree)" : "x"))
end

########################
# Queries about a term #
########################

"""
Check if two terms are equal. (Should this check if the coefficient types are equal?)
"""

==(t1::Term,t2::Term)::Bool = (t1.coeff == t2.coeff) && (t1.degree == t2.degree)

"""
Check if a term is 0.
"""
iszero(t::Term)::Bool = iszero(t.coeff)

"""
Compare two terms.
"""
isless(t1::Term,t2::Term)::Bool =  t1.degree == t2.degree ? (t1.coeff < t2.coeff) : (t1.degree < t2.degree)  

"""
Evaluate a term at a point x.
"""
evaluate(t::Term, x::T) where T <: Number = t.coeff * x^t.degree

##########################
# Operations with a term #
##########################

"""
Add two terms of the same degree.
"""
function +(t1::Term,t2::Term)::Term
    @assert t1.degree == t2.degree
    Term(t1.coeff + t2.coeff, t1.degree)
end

"""
Negate a term.
"""
-(t::Term,) = Term(-t.coeff,t.degree)  

"""
Subtract two terms with the same degree.
"""
-(t1::Term, t2::Term)::Term = t1 + (-t2) 

"""
Multiply two terms.
"""
*(t1::Term, t2::Term)::Term = Term(t1.coeff * t2.coeff, t1.degree + t2.degree)


"""
Compute the symmetric mod of a term with an integer.
"""
mod(t::Term, p::Integer) = Term(mod(t.coeff,p), t.degree)

"""
Compute the derivative of a term.
"""
derivative(t::Term) = Term(t.coeff*t.degree,max(t.degree-1,0))

"""
Divide two terms. Returns a function of an integer.
"""
function ÷(t1::Term,t2::Term) #\div + [TAB]
    @assert t1.degree ≥ t2.degree
    f(p::Integer)::Term = Term(mod((t1.coeff * int_inverse_mod(t2.coeff, p)), p), t1.degree - t2.degree)
end

"""
Integer divide a term by an integer.
"""
÷(t::Term, n::Integer) = t ÷ Term(n,0)
