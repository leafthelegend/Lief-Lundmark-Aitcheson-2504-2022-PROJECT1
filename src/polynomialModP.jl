struct PolynomialModP
    polynomial::Polynomial
    p::Integer

    function PolynomialModP(polynomial::Polynomial,p::Integer;reducing=true)
        if reducing
            return new(mod(polynomial,p),p)
        else
            return new(polynomial,p)
        end
    end
end

#inplace polynomial mod
function mod!(f::Polynomial, p::Integer)::Polynomial
    for i in 1:length(f.terms)
        f.terms[i] = mod(f.terms[i], p)
    end
    return trim!(f)
end

function reduce!(polynomial::PolynomialModP)::PolynomialModP
    mod!(polynomial.polynomial, polynomial.p)
    return polynomial
end

#display
function show(io::IO,p::PolynomialModP)
    show(io,p.polynomial)
    print(io," (mod $(p.p))")
end

#derivative
function derivative(p::PolynomialModP)
    return PolynomialModP(derivative(p.polynomial),p.p)
end

#arithmetic methods

#addition

#add a polynomialModP and a term:
+(p::PolynomialModP, t::Term)::PolynomialModP = PolynomialModP(p.polynomial+t,p.p)

+(t::Term,p::PolynomialModP) = p + t

#add a polynomialModP and a polynomial:
+(pm::PolynomialModP, p::Polynomial)::PolynomialModP = PolynomialModP(pm.polynomial+p,pm.p)

+(p::Polynomial,pm::PolynomialModP) = pm + p

#add a polynomialModP and an integer:
+(p::PolynomialModP, n::Integer)::PolynomialModP = PolynomialModP(p.polynomial+n,p.p)

+(n::Integer,p::PolynomialModP) = p + n

#add two polynomialModPs:
function +(p1::PolynomialModP,p2::PolynomialModP)::PolynomialModP
    @assert p1.p == p2.p
    return PolynomialModP(p1.polynomial+p2.polynomial,p1.p)
end

#take the negative of a polynomialModP: #note - why does this not work with a return type annotation???
-(p::PolynomialModP) = PolynomialModP(-p.polynomial,p.p)

#Subtraction
-(p1::PolynomialModP,p2::PolynomialModP)::PolynomialModP = p1 + (-p2)

-(p1::PolynomialModP,p::Polynomial)::PolynomialModP = p1 + (-p)

-(p::Polynomial,p1::PolynomialModP)::PolynomialModP = p + (-p1)

-(p1::PolynomialModP,t::Term)::PolynomialModP = p1 + (-t)

-(t::Term,p1::PolynomialModP)::PolynomialModP = t + (-p1)

-(p1::PolynomialModP,n::Integer)::PolynomialModP = p1 + (-n)

-(n::Integer,p1::PolynomialModP)::PolynomialModP = n + (-p1)

#multiplication

#multiply a polynomialModP and a term:
*(p::PolynomialModP, t::Term)::PolynomialModP = PolynomialModP(p.polynomial*t,p.p)

*(t::Term,p::PolynomialModP) = p * t

#multiply a polynomialModP and a polynomial:
*(pm::PolynomialModP, p::Polynomial)::PolynomialModP = PolynomialModP(pm.polynomial*p,pm.p)

*(p::Polynomial,pm::PolynomialModP) = pm * p

#multiply a polynomialModP and an integer:
*(p::PolynomialModP, n::Integer)::PolynomialModP = PolynomialModP(p.polynomial*n,p.p)

*(n::Integer,p::PolynomialModP) = p * n

#multiply two polynomialModPs:
function *(p1::PolynomialModP,p2::PolynomialModP)::PolynomialModP
    @assert p1.p == p2.p
    return PolynomialModP(p1.polynomial*p2.polynomial,p1.p)
end

#exponentiation
# function ^(p::PolynomialModP,n::Integer)::PolynomialModP
#     n < 0 && error("No negative power")
#     out = p
#     for _ in 1:n-1
#         out *=p
#     end
#     return out
# end


#divide a polynomialModP by an integer:
÷(p::PolynomialModP, n::Integer)::PolynomialModP = PolynomialModP((p.polynomial÷n)(p.p),p.p)

#divide two polynomialModPs:

function divide(num::PolynomialModP, den::PolynomialModP)
    @assert num.p == den.p
    return PolynomialModP.(divide(num.polynomial,den.polynomial)(num.p),num.p)
end

function ÷(num::PolynomialModP, den::PolynomialModP)::PolynomialModP
    @assert num.p == den.p
    return PolynomialModP((num.polynomial÷den.polynomial)(num.p),num.p)
end

function rem(num::PolynomialModP, den::PolynomialModP)::PolynomialModP
    @assert num.p == den.p
    return PolynomialModP(rem(num.polynomial,den.polynomial)(num.p),num.p)
end

#gcd of two polynomialModPs:
function gcd(a::PolynomialModP, b::PolynomialModP)::PolynomialModP
    @assert a.p == b.p
    return PolynomialModP(gcd(a.polynomial,b.polynomial,a.p),a.p)
end

# prim_part(a::PolynomialModP) = PolynomialModP(prim_part(a.polynomial)(a.p),a.p)
prim_part(a::PolynomialModP) = a ÷ content(a.polynomial)

# #factorise a polynomialModP:
# function factor(p::PolynomialModP)::Vector{Tuple{PolynomialModP,Integer}}
#     factors = factor(p.polynomial,p.p)
#     return map(f->(PolynomialModP(f[1],p.p),f[2]),factors)
# end

# function expand_factorization(factors::Vector{<:Tuple{PolynomialModP,Integer}})
#     return PolynomialModP(expand_factorization(map(f->(f[1].polynomial,f[2]),factors)),factors[1][1].p)
# end

#boolean methods:
==(p1::PolynomialModP,p2::PolynomialModP) = p1.p == p2.p && p1.polynomial == p2.polynomial 

iszero(p::PolynomialModP) = iszero(p.polynomial)

==(p::PolynomialModP, n::Integer) = iszero(p - n)

#creation
zero(p::PolynomialModP) = PolynomialModP(zero(p.polynomial),p.p)
one(p::PolynomialModP) = PolynomialModP(one(p.polynomial),p.p)
x_poly(p::PolynomialModP) = PolynomialModP(x_poly(typeof(p.polynomial)),p.p)
x_poly_modp(polyType::Type{<:Polynomial},p::Int) = PolynomialModP(x_poly(polyType),p)
copyTerms(p::PolynomialModP) = PolynomialModP(copyTerms(p.polynomial),p.p,reducing=false)
rand_modp(polyType::Type{<:Polynomial},p::Int;kwargs...) = PolynomialModP(rand(polyType;kwargs...),p)

#queries
degree(p::PolynomialModP) = degree(p.polynomial)
leading(p::PolynomialModP) = leading(p.polynomial)
