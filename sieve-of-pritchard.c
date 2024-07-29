#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>

void sieve_of_pritchard(mpz_t limit) {
    mpz_t n, m, a, b, c, i, j, p, q, x, y, z;
    mpz_inits(n, m, a, b, c, i, j, p, q, x, y, z, NULL);
    mpz_set_ui(n, 1);
    mpz_set_ui(m, 1);
    mpz_set_ui(a, 1);
    mpz_set_ui(b, 1);
    mpz_set_ui(c, 1);
    
    mpz_t *primes = malloc(mpz_get_ui(limit) * sizeof(mpz_t));
    int *is_prime = malloc(mpz_get_ui(limit) * sizeof(int));
    for (mpz_set_ui(i, 2); mpz_cmp(i, limit) <= 0; mpz_add_ui(i, i, 1)) {
        mpz_init(primes[mpz_get_ui(i)]);
        is_prime[mpz_get_ui(i)] = 1;
    }

    for (mpz_set_ui(i, 2); mpz_cmp(i, limit) <= 0; mpz_add_ui(i, i, 1)) {
        if (is_prime[mpz_get_ui(i)]) {
            mpz_set(primes[mpz_get_ui(m)], i);
            mpz_add_ui(m, m, 1);
            for (mpz_mul(j, i, i); mpz_cmp(j, limit) <= 0; mpz_add(j, j, i)) {
                is_prime[mpz_get_ui(j)] = 0;
            }
        }
    }

    for (mpz_set_ui(i, 0); mpz_cmp(i, m) < 0; mpz_add_ui(i, i, 1)) {
        gmp_printf("%Zd\n", primes[mpz_get_ui(i)]);
    }

    for (mpz_set_ui(i, 2); mpz_cmp(i, limit) <= 0; mpz_add_ui(i, i, 1)) {
        mpz_clear(primes[mpz_get_ui(i)]);
    }

    free(primes);
    free(is_prime);
    mpz_clears(n, m, a, b, c, i, j, p, q, x, y, z, NULL);
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: %s <limit>\n", argv[0]);
        return 1;
    }

    mpz_t limit;
    mpz_init(limit);
    mpz_set_str(limit, argv[1], 10);
    sieve_of_pritchard(limit);
    mpz_clear(limit);

    return 0;
}
