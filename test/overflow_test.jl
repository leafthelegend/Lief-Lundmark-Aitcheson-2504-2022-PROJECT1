function test_overflow(T::Type{<:Polynomial};n::Int=128) #n=128 will overflow for all Julia's primitive integer types
    p = 2*one(T)
    res = one(T)
    for _ in 1:n
        newRes = res*p
        @assert leading(newRes).coeff > leading(res).coeff
        res = newRes
    end
    println("test_overflow - PASSED")
end