function nth_prime(n)
   primes = Int[]
   candidate = 2
   while length(primes) < n
       if all(candidate % prime != 0 for prime in primes)
           push!(primes, candidate)
       end
       candidate += 1
   end
   return primes[end]
end

println(nth_prime(10)) # prints 29
