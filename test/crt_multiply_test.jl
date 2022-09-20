power = 64
mean_deg = 5000.

function slow_prod_test(;N::Int = 10, N_prods::Int = 0, seed::Int = 0, verbose::Bool = true)
    Random.seed!(seed)
    verbose && println("Multiplying pairs:")
    for _ in 1:N
        p1 = rand(PolynomialSparse{BigInt},max_coeff=big(2)^power,mean_degree=mean_deg)
        p2 = rand(PolynomialSparse{BigInt},max_coeff=big(2)^power,mean_degree=mean_deg)
        prod = p1*p2
        @assert leading(prod) == leading(p1)*leading(p2)
        verbose && print(".")
    end
    verbose && print("\nMultiplying chains:")
    for _ in 1:1
        verbose && print("\n|")
        p_base = PolynomialSparse{BigInt}(Term{BigInt}(1,0))
        for _ in 1:N_prods
            p = rand(PolynomialSparse{BigInt},max_coeff=big(2)^power,mean_degree=mean_deg)
            prod = p_base*p
            #we must also handle the case where the leading term is zero mod p, so the first part of this assertion will fail
            @assert leading(prod) == leading(p_base)*leading(p)
            p_base = prod
            verbose && print(".")
        end
    end
    verbose && println("\nslow_prod_test - PASSED")
    nothing
end

function crt_prod_test(;N::Int = 10, N_prods::Int = 0, seed::Int = 0, verbose::Bool = true)
    Random.seed!(seed)
    verbose && println("Multiplying pairs:")
    for _ in 1:N
        p1 = rand(PolynomialSparse{BigInt},max_coeff=big(2)^power,mean_degree=mean_deg)
        p2 = rand(PolynomialSparse{BigInt},max_coeff=big(2)^power,mean_degree=mean_deg)
        prod = crt_multiply(p1,p2)
        @assert leading(prod) == leading(p1)*leading(p2)
        verbose && print(".")
    end
    verbose && print("\nMultiplying chains:")
    for _ in 1:1
        verbose && print("\n|")
        p_base = PolynomialSparse{BigInt}(Term{BigInt}(1,0))
        for _ in 1:N_prods
            p = rand(PolynomialSparse{BigInt},max_coeff=big(2)^power,mean_degree=mean_deg)
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