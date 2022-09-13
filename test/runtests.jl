#############################################################################
#############################################################################
#
# A script that runs all unit tests in the project.
#                                                                               
#############################################################################
#############################################################################

using Pkg
Pkg.activate(".")

include("../poly_factorization_project.jl")

polyType = PolynomialSparse{BigInt} #test sparse or dense

####
# Execute unit tests for integers
###
include("integers_test.jl")
test_euclid_ints()
test_ext_euclid_ints()

####
# Execute unit tests for polynomials
####
include("polynomials_test.jl")
# prod_test_poly(polyType)
# prod_derivative_test_poly(polyType)
ext_euclid_test_poly(polyType)
division_test_poly(polyType)

####
# Execute unit tests for polynomial factorization
####
include("factorization_test.jl")
factor_test_poly(polyType)