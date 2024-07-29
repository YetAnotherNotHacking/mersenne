# Mersenne Cache: 
Mersenne with cache to avoid filling memory, which makes the gnu arithmetic library even more useful.
This is the best script in this repo for launching it and leaving for a while, unless of course you 
are on datacenter hardware, then use the intel script.

# Mersenne:
standard mersenne prime search, without cache to disk. Works with large amounts of memory on machines
that have under 32 threads (16 cores if hyperthreading is enabled)

# Prime:
slow standard expanding prime search with support for setting and initial length of prime to start at 
or something like that, it is old script, and not optimized in any way shape or form.



# Mersenne-intel:
a script for the intel cloud developer platform, it is able to find very very large 
primes on this cpu, the high seirra that they allow people to test on in the cloud for free. since there 
is such a large number of cores, I made my own "super cores"

I say cores a lot here, but I mean hyperthreads (linux sees them as cores)

the super cores are working as so: 16 cores turns into one super core, do the math and there are 24 super cores, and the program is just mersenne adapted ot be able to group cores together and not have them overwhelm the program

the more cores there are, the more work it is just to use the cores effeciently, but this is able to run this on special cores that do not overwhelm, and are made to work together in a more effecient way.

# Seive of Eratoshtenes
https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
don't got no reason to yap about it here
