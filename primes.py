def nth_prime(n):
   primes = [2, 3]
   i = 3
   while len(primes) < n:
       if all(i % prime != 0 for prime in primes):
           primes.append(i)
       i += 1
   return primes[-1]


print(nth_prime(10))
