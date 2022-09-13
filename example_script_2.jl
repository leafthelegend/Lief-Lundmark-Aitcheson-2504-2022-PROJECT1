using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

x = x_poly(PolynomialSparse)
p1 = 2x^3 + 4x^2 - 3x
p2 = 2x^4 - 4x^2 - 3x + 3

println("We can create polynomials and perform arithmetic with them. Here we create two polynomials, p1 and p2:\n")

println("p1 = $p1")
println("p2 = $p2\n")

println("We can compute their sum and product, or integer powers:\n")

println("p1+p2 = $(p1+p2)")
println("p1 × p2 = $(p1*p2)")
println("p2^4 = $(p2^4)\n")

println("We can also compute derivatives:\n")

println("d/dx p2 = $(derivative(p2))")
println("d/dx (p1*p2) = $(derivative(p1*p2))\n")

println("If we work with polynomials over prime residue classes, we can perform polynomial division and factorisation.")

prime = 17
p = mod((7x^3 + 2x^2 + 8x + 1)*(x^2+x+1),prime)

println("Consider the polynomial $(p) (mod $prime).")
factorisation = factor(p,prime)
println("We can factorise this as:")
print.("($(f[1]))$(f[2]>1 ? "^$(f[1])" : "")" for f in factorisation)
println()
println("We can also go the opposite direction, and expand the factorised form.")

pr = mod(expand_factorization(factorisation),prime)
println("Expanded polynomial (mod $prime) = ", pr, "\n")
println("We can use the extended euclidean algorithm to compute the GCD and Bezout coefficients of a pair of polynomials.")

p3 = 99x^3+45x^2+77x+82
p4 = 98x^4+39x^2+97x+15
println("Consider the following polynomials, p3 and p4 (mod 101):")
println("p3 = $p3")
println("p4 = $p4")

g,s,t = extended_euclid_alg(p3,p4,101)

println("We can compute gcd(p3,p4) to be $g (mod 101).")
println("The Bezout coefficients are then polynomials s and t such that s*p3 + t*p4 = $g.")
println("We compute these polynomials to be:")
println("s = $s")
println("t = $t")
