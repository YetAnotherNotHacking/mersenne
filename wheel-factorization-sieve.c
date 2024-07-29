#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gmp.h>

#define WHEEL_SIZE 48
#define WHEEL_PRIMES 5

const unsigned long wheel[WHEEL_SIZE] = {
    1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 
    67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 121, 127, 131, 137, 
    139, 143, 149, 151, 157, 163, 167, 169, 173, 179, 181, 187, 191, 193, 197, 199
};

bool is_prime(mpz_t n, mpz_t *primes, int prime_count) {
    if (mpz_cmp_ui(n, 2) < 0) return false;
    
    mpz_t sqrt_n, temp;
    mpz_init(sqrt_n);
    mpz_init(temp);
    
    mpz_sqrt(sqrt_n, n);
    
    for (int i = 0; i < prime_count; i++) {
        if (mpz_cmp(primes[i], sqrt_n) > 0) break;
        mpz_mod(temp, n, primes[i]);
        if (mpz_cmp_ui(temp, 0) == 0) {
            mpz_clear(sqrt_n);
            mpz_clear(temp);
            return false;
        }
    }
    
    mpz_clear(sqrt_n);
    mpz_clear(temp);
    return true;
}

void find_primes(mpz_t limit) {
    mpz_t *primes = malloc(sizeof(mpz_t) * mpz_get_ui(limit));
    int prime_count = 0;

    // Add first few primes manually
    unsigned long first_primes[WHEEL_PRIMES] = {2, 3, 5, 7, 11};
    for (int i = 0; i < WHEEL_PRIMES; i++) {
        mpz_t temp;
        mpz_init_set_ui(temp, first_primes[i]);
        if (mpz_cmp(temp, limit) <= 0) {
            mpz_init_set(primes[prime_count], temp);
            gmp_printf("%Zd ", primes[prime_count]);
            prime_count++;
        }
        mpz_clear(temp);
    }

    mpz_t n;
    mpz_init_set_ui(n, 13);

    while (mpz_cmp(n, limit) <= 0) {
        if (is_prime(n, primes, prime_count)) {
            mpz_init_set(primes[prime_count], n);
            gmp_printf("%Zd ", primes[prime_count]);
            prime_count++;
        }
        mpz_add_ui(n, n, 2);
    }

    for (int i = 0; i < prime_count; i++) {
        mpz_clear(primes[i]);
    }
    free(primes);
    mpz_clear(n);
    printf("\n");
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <upper_limit>\n", argv[0]);
        return 1;
    }

    mpz_t limit;
    mpz_init(limit);

    if (mpz_set_str(limit, argv[1], 10) != 0) {
        fprintf(stderr, "Error: Invalid number format. Please provide a positive integer.\n");
        mpz_clear(limit);
        return 1;
    }

    if (mpz_cmp_ui(limit, 0) <= 0) {
        fprintf(stderr, "Error: Please provide a positive integer as the upper limit.\n");
        mpz_clear(limit);
        return 1;
    }

    gmp_printf("Prime numbers up to %Zd are:\n", limit);
    find_primes(limit);

    mpz_clear(limit);
    return 0;
}
