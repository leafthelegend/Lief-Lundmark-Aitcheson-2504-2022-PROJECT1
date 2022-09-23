function slow_mult(p1::Polynomial, p2::Polynomial)::Polynomial
    p_out = typeof(p1)()
    for t in p1
        new_summand = (t * p2)
        p_out = add!(p_out, new_summand)
    end
    return p_out
end

function slow_prod_test(;N::Int = 10, N_prods::Int = 0, seed::Int = 0, verbose::Bool = true, deg::Int = 5000,power::Int=512)
    Random.seed!(seed)
    verbose && println("Multiplying pairs:")
    for _ in 1:N
        p1 = rand(PolynomialSparse{BigInt},max_coeff=big(2)^power,degree=deg)
        p2 = rand(PolynomialSparse{BigInt},max_coeff=big(2)^power,degree=deg)
        prod = slow_mult(p1,p2)
        @assert leading(prod) == leading(p1)*leading(p2)
        verbose && print(".")
    end
    verbose && print("\nMultiplying chains:")
    for _ in 1:N
        verbose && print("\n|")
        p_base = PolynomialSparse{BigInt}(Term{BigInt}(1,0))
        for _ in 1:N
            p = rand(PolynomialSparse{BigInt},max_coeff=big(2)^power,degree=deg)
            prod = slow_mult(p_base,p)
            #we must also handle the case where the leading term is zero mod p, so the first part of this assertion will fail
            @assert leading(prod) == leading(p_base)*leading(p)
            p_base = prod
            verbose && print(".")
        end
    end
    verbose && println("\nslow_prod_test - PASSED")
    nothing
end

function crt_prod_test(;N::Int = 10, N_prods::Int = 0, seed::Int = 0, verbose::Bool = true, deg::Int = 5000, power::Int=512)
    Random.seed!(seed)
    verbose && println("Multiplying pairs:")
    for _ in 1:N
        p1 = rand(PolynomialSparse{BigInt},max_coeff=big(2)^power,degree=deg)
        p2 = rand(PolynomialSparse{BigInt},max_coeff=big(2)^power,degree=deg)
        prod = crt_multiply(p1,p2)
        @assert leading(prod) == leading(p1)*leading(p2)
        verbose && print(".")
    end
    verbose && print("\nMultiplying chains:")
    for _ in 1:N
        verbose && print("\n|")
        p_base = PolynomialSparse{BigInt}(Term{BigInt}(1,0))
        for _ in 1:N_prods
            p = rand(PolynomialSparse{BigInt},max_coeff=big(2)^power,degree=deg)
            prod = crt_multiply(p_base,p)
            #we must also handle the case where the leading term is zero mod p, so the first part of this assertion will fail
            @assert leading(prod) == leading(p_base)*leading(p)
            p_base = prod
            verbose && print(".")
        end
    end
    verbose && println("\ncrt_prod_test - PASSED")
    nothing
end

function verify_crt(;N::Int = 100, seed::Int = 0, verbose::Bool = true, mean_deg::Float64 = 100., power::Int=20)
    Random.seed!(seed)
    verbose && println("Verifying CRT:")
    for _ in 1:N
        p1 = rand(PolynomialSparse{BigInt},max_coeff=big(2)^power,mean_degree=mean_deg)
        p2 = rand(PolynomialSparse{BigInt},max_coeff=big(2)^power,mean_degree=mean_deg)
        prod = crt_multiply(p1,p2)
        @assert prod == slow_mult(p1,p2)
        verbose && print(".")
    end
    verbose && println("\nverify_crt - PASSED")
    nothing
end