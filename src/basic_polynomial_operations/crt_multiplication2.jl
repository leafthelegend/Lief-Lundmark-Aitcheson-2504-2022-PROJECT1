function CRT(a::Polynomial,b::Polynomial,p1::Integer,p2::Integer)::Polynomial
    T = typeof(a)
    v1 = copyTerms(a)
    v2 = b-v1
    v2*=int_inverse_mod(p1,p2)
    v2 = mod!(v2,p2)
    p=add!(v1,v2*p1)
    @assert mod(p,p1)==a && mod(p,p2)==b
    return p
end

#convert to and from big int and Int64 polynomials

smallTerm(t::Term{BigInt})::Term{Int64} = Term{Int64}(Int64(t.coeff),t.degree)
smallPoly(p::PolynomialSparse{BigInt})::PolynomialSparse{Int64} = PolynomialSparse{Int64}(smallTerm.(p.terms))

bigTerm(t::Term{Int64})::Term{BigInt} = Term{BigInt}(BigInt(t.coeff),t.degree)
bigPoly(p::PolynomialSparse{Int64})::PolynomialSparse{BigInt} = PolynomialSparse{BigInt}(bigTerm.(p.terms))

smallPoly(p::PolynomialModP)::PolynomialModP = PolynomialModP(smallPoly(p.polynomial),p.p)
bigPoly(p::PolynomialModP)::PolynomialModP = PolynomialModP(bigPoly(p.polynomial),p.p)

using Primes

function crt_multiply(p::PolynomialSparse{BigInt},q::PolynomialSparse{BigInt})
    # println("p=$p")
    # println("q=$q")
    height = 2*max(abs.(coeffs(p))...)*max(abs.(coeffs(q))...)*min(length(p),length(q))

    # primes = [Primes.prevprime(typemax(Int))]
    # M=big(primes[1])
    # while M <= height
    #     prime::Int64 = Primes.prevprime(primes[1]-1)
    #     M *= prime
    #     pushfirst!(primes,prime)
    # end
    prime::Int = Primes.prevprime(2^31)
    M=big(prime)
    c = bigPoly((smallPoly(PolynomialModP(p,prime))*smallPoly(PolynomialModP(q,prime))).polynomial)
    while M <= height
        prime = Primes.prevprime(prime-1)
        @assert prime>3
        c_new = smallPoly(PolynomialModP(p,prime))*smallPoly(PolynomialModP(q,prime))
        c = CRT(c,bigPoly(c_new.polynomial),M,prime)
        M *= prime
    end

    return smod(c,M)
end