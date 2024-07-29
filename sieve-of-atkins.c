#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

void sieve_of_atkin(mpz_t limit) {
    mpz_t x2, y2, n, limit_sqrt, remainder, temp, k;
    mpz_init(x2);
    mpz_init(y2);
    mpz_init(n);
    mpz_init(limit_sqrt);
    mpz_init(remainder);
    mpz_init(temp);
    mpz_init(k);

    int limit_int = mpz_get_ui(limit);
    char* sieve = calloc(limit_int + 1, sizeof(char));

    mpz_sqrt(limit_sqrt, limit);
    for (mpz_set_ui(n, 1); mpz_cmp(n, limit_sqrt) <= 0; mpz_add_ui(n, n, 1)) {
        mpz_mul(x2, n, n);

        for (mpz_set_ui(y2, 1); mpz_cmp(y2, limit_sqrt) <= 0; mpz_add_ui(y2, y2, 1)) {
            mpz_mul(y2, y2, y2);

            mpz_set(temp, x2);
            mpz_add(temp, temp, y2);
            if (mpz_cmp(temp, limit) <= 0) {
                mpz_mod_ui(remainder, temp, 12);
                int r = mpz_get_ui(remainder);
                if (r == 1 || r == 5) {
                    sieve[mpz_get_ui(temp)] ^= 1;
                }
            }

            mpz_set(temp, x2);
            mpz_mul_ui(temp, temp, 3);
            mpz_sub(temp, temp, y2);
            if (mpz_cmp_ui(temp, 0) > 0 && mpz_cmp(temp, limit) <= 0) {
                mpz_mod_ui(remainder, temp, 12);
                int r = mpz_get_ui(remainder);
                if (r == 7) {
                    sieve[mpz_get_ui(temp)] ^= 1;
                }
            }

            mpz_set(temp, x2);
            mpz_mul_ui(temp, temp, 3);
            mpz_add(temp, temp, y2);
            if (mpz_cmp(temp, limit) <= 0) {
                mpz_mod_ui(remainder, temp, 12);
                int r = mpz_get_ui(remainder);
                if (r == 11) {
                    sieve[mpz_get_ui(temp)] ^= 1;
                }
            }
        }
    }

    for (mpz_set_ui(n, 5); mpz_cmp(n, limit_sqrt) <= 0; mpz_add_ui(n, n, 1)) {
        if (sieve[mpz_get_ui(n)]) {
            mpz_mul(temp, n, n);
            for (mpz_set(k, temp); mpz_cmp(k, limit) <= 0; mpz_add(k, k, temp)) {
                sieve[mpz_get_ui(k)] = 0;
            }
        }
    }

    if (limit_int >= 2) printf("2 ");
    if (limit_int >= 3) printf("3 ");
    for (mpz_set_ui(n, 5); mpz_cmp(n, limit) <= 0; mpz_add_ui(n, n, 1)) {
        if (sieve[mpz_get_ui(n)]) {
            gmp_printf("%Zd ", n);
        }
    }
    printf("\n");

    mpz_clear(x2);
    mpz_clear(y2);
    mpz_clear(n);
    mpz_clear(limit_sqrt);
    mpz_clear(remainder);
    mpz_clear(temp);
    mpz_clear(k);
    free(sieve);
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <limit>\n", argv[0]);
        return 1;
    }

    mpz_t limit;
    mpz_init_set_str(limit, argv[1], 10);

    sieve_of_atkin(limit);

    mpz_clear(limit);
    return 0;
}
