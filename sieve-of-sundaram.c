// used chatgpt4o to get a base for a threaded version, but this was not easy to make threaded. I changed it myself to not be threaded and got it working

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

void sieve_of_sundaram(int n, int *marked) {
    int half = (n - 1) / 2;
    for (int i = 1; i <= half; ++i) {
        for (int j = i; (i + j + 2 * i * j) <= half; ++j) {
            marked[i + j + 2 * i * j] = 1;
        }
    }
}

void print_primes(int n, int *marked) {
    mpz_t prime;
    mpz_init(prime);

    if (n > 2) {
        printf("2\n");
    }

    for (int i = 1; i <= (n - 1) / 2; ++i) {
        if (!marked[i]) {
            mpz_set_ui(prime, 2 * i + 1);
            gmp_printf("%Zd\n", prime);
        }
    }

    mpz_clear(prime);
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <upper_bound>\n", argv[0]);
        return EXIT_FAILURE;
    }

    int n = atoi(argv[1]);
    if (n < 2) {
        fprintf(stderr, "Upper bound must be at least 2\n");
        return EXIT_FAILURE;
    }

    int *marked = calloc((n - 1) / 2 + 1, sizeof(int));
    if (marked == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return EXIT_FAILURE;
    }

    sieve_of_sundaram(n, marked);
    print_primes(n, marked);

    free(marked);
    return EXIT_SUCCESS;
}
