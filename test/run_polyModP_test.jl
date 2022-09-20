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

polyType = PolynomialSparse{BigInt} #which polynomial type to test


####
# Execute unit tests for polynomials
####
include("polyModP_test.jl")
prod_test_polyModP(polyType)
division_test_polyModP(polyType)
factor_test_polyModP(polyType)