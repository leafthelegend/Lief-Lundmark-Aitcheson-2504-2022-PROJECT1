#############################################################################
#############################################################################
#
# This is the main project file for polynomial factorization
#                                                                               
#############################################################################
#############################################################################
import Pkg;


using Distributions, StatsBase, Random

import Base: %
import Base: push!, pop!, iszero, show, isless, map, map!, iterate, length, last
import Base: +, -, *, mod, %, ÷, ==, ^, rand, rem, zero, one

include("src/general_alg.jl")
include("src/term.jl")
include("src/polynomial.jl")
include("src/polynomialDense.jl")
include("src/polynomialSparse.jl")
include("src/polynomialModP.jl")
include("src/basic_polynomial_operations/polynomial_addition.jl")
include("src/basic_polynomial_operations/polynomial_multiplication.jl")
include("src/basic_polynomial_operations/powers.jl")
include("src/basic_polynomial_operations/polynomial_division.jl")
include("src/basic_polynomial_operations/polynomial_gcd.jl")
include("src/polynomial_factorization/factorModP.jl")
include("src/basic_polynomial_operations/crt_multiplication2.jl")

nothing