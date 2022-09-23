struct PolynomialSparse{T<:Integer} <: Polynomial

    #A vector of non-zero terms, in descending order
    terms::Vector{Term{T}}

    #Inner constructor of 0 polynomial
    PolynomialSparse{T}() where {T<:Integer} = new{T}([])
    PolynomialSparse() = new{Int}([])

    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialSparse{T}(vt::Vector{Term{T}};sorting=true,filtering=true) where {T<:Integer}

        #Remove zero terms
        if filtering
            vt = filter((t)->!iszero(t), vt)
        end

        #sort the terms in descending order of degree
        if sorting
            vt = sort(vt,by=t->t.degree,rev=true)
        end
        return new{T}(vt)
    end
end

PolynomialSparse(vt::Vector{Term{T}}) where {T<:Integer} = PolynomialSparse{T}(vt)

#alias
const PolynomialSparseBI = PolynomialSparse{BigInt}

#implement parametric methods to handle the case where a coefficient type is supplied

cyclotonic_polynomial(polyType::Type{<:PolynomialSparse{T}},p::Int) where {T<:Integer} = polyType([Term{T}(1,p), Term{T}(-1,0)])

linear_monic_polynomial(polyType::Type{<:PolynomialSparse{T}},n::Int) where {T<:Integer} = polyType([Term{T}(1,1), Term{T}(-n,0)])

x_poly(polyType::Type{<:PolynomialSparse{T}}) where {T<:Integer} = polyType(Term{T}(1,1))

one(polyType::Type{<:PolynomialSparse{T}}) where {T<:Integer} = polyType(one(Term{T}))

function rand(polyType::Type{PolynomialSparse{T}} ; 
    degree::Int = -1, 
    terms::Int = -1, 
    max_coeff::T = T(100), 
    mean_degree::Float64 = 5.0,
    prob_term::Float64  = 0.7,
    monic = false,
    condition = (p)->true) where {T<:Integer}

while true 
_degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
_terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
coeffs = rand(-max_coeff:max_coeff,_terms+1)
monic && (coeffs[end] = 1)
p = polyType( [Term{T}(coeffs[i],degrees[i]) for i in 1:length(degrees)] )
condition(p) && return p
end
end

"""
This function maintains the invariant of the PolynomialSparse type so that there are no zero terms.
"""
function trim!(p::PolynomialSparse)::PolynomialSparse
    filter!((t)->!iszero(t), p.terms)
    return p
end

"""
Check if the polynomial is zero.
"""
iszero(p::PolynomialSparse)::Bool = isempty(p.terms)

"""
Construct a sparse polynomial with a single term.
"""
PolynomialSparse{T}(t::Term) where {T<:Integer}= PolynomialSparse{T}([t])
PolynomialSparse(t::Term) = PolynomialSparse([t])

###########
# Display #
###########

"""
Show a sparse polynomial.
"""
function show(io::IO, p::PolynomialSparse)
    if iszero(p)
        print(io,"0")
    else
        n = length(p.terms)
        for (i,t) in enumerate(p.terms)
                print(io, i == 1 ? t : "$(t.coeff>0 ? "+" : "")$(t)")
        end
    end
end

##############################
# Queries about a sparse polynomial #
##############################

"""
The leading term of the polynomial.
"""
leading(p::PolynomialSparse)::Term = isempty(p.terms) ? zero(Term) : first(p.terms) 

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the polynomial.
"""
#Note that ideally this would throw and error if pushing another term of degree that is already in the polynomial
function push!(p::PolynomialSparse, t::Term)
    if isempty(p.terms)
        push!(p.terms,t)
        return p
    end
    i = findfirst(term -> term.degree <= t.degree, p.terms)
    if isnothing(i)
        push!(p.terms,t)
    elseif t.degree == p.terms[i].degree
        p.terms[i] = t
    else
        insert!(p.terms,i,t)
    end
    return p     
end

"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.
"""
function pop!(p::PolynomialSparse)::Term 
    if iszero(p)
        return zero(Term)
    end
    return popfirst!(p.terms)
end
