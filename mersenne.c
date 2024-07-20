#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <gmp.h>
#include <signal.h>
#include <unistd.h>
#include <time.h>
#include <getopt.h>

#define MAX_THREADS 64
#define MILLER_RABIN_ITERATIONS 40
#define UPDATE_INTERVAL 1000 // Update status every 1000 checks

// ANSI color codes
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

typedef struct {
    mpz_t start;
    mpz_t step;
    int thread_id;
} thread_data_t;

volatile sig_atomic_t keep_running = 1;
pthread_mutex_t prime_mutex = PTHREAD_MUTEX_INITIALIZER;
mpz_t current_prime;
unsigned long long current_n;
int num_threads;
unsigned long long primes_checked = 0;
time_t start_time;

// Function prototype
void print_status(void);

void handle_sigint(int sig) {
    keep_running = 0;
}

int miller_rabin(mpz_t n, int iterations) {
    if (mpz_cmp_ui(n, 2) < 0)
        return 0;
    if (mpz_cmp_ui(n, 2) == 0)
        return 1;
    if (mpz_even_p(n))
        return 0;

    mpz_t a, y, j, r, n_minus_one;
    mpz_inits(a, y, j, r, n_minus_one, NULL);

    unsigned long int s = 0;
    mpz_sub_ui(n_minus_one, n, 1);
    mpz_set(r, n_minus_one);

    while (mpz_even_p(r)) {
        mpz_fdiv_q_2exp(r, r, 1);
        s++;
    }

    gmp_randstate_t rnd;
    gmp_randinit_default(rnd);
    gmp_randseed_ui(rnd, time(NULL));

    int is_prime = 1;
    for (int i = 0; i < iterations && is_prime; i++) {
        mpz_urandomm(a, rnd, n_minus_one);
        mpz_add_ui(a, a, 1);
        mpz_powm(y, a, r, n);

        if (mpz_cmp_ui(y, 1) != 0 && mpz_cmp(y, n_minus_one) != 0) {
            for (unsigned long int j = 1; j < s && mpz_cmp(y, n_minus_one) != 0; j++) {
                mpz_powm_ui(y, y, 2, n);
                if (mpz_cmp_ui(y, 1) == 0) {
                    is_prime = 0;
                    break;
                }
            }
            if (mpz_cmp(y, n_minus_one) != 0) {
                is_prime = 0;
            }
        }
    }

    mpz_clears(a, y, j, r, n_minus_one, NULL);
    gmp_randclear(rnd);
    return is_prime;
}

void* find_mersenne_primes(void* arg) {
    thread_data_t* data = (thread_data_t*)arg;
    mpz_t candidate, mersenne;
    mpz_inits(candidate, mersenne, NULL);
    mpz_set(candidate, data->start);

    while (keep_running) {
        mpz_ui_pow_ui(mersenne, 2, mpz_get_ui(candidate));
        mpz_sub_ui(mersenne, mersenne, 1);

        if (miller_rabin(mersenne, MILLER_RABIN_ITERATIONS)) {
            pthread_mutex_lock(&prime_mutex);
            if (mpz_cmp(candidate, current_prime) > 0) {
                mpz_set(current_prime, candidate);
                current_n = mpz_get_ui(candidate);
                printf(ANSI_COLOR_GREEN "\nFound Mersenne prime: 2^%llu - 1\n" ANSI_COLOR_RESET, current_n);
                mpz_out_str(stdout, 10, mersenne);
                printf("\n");
            }
            pthread_mutex_unlock(&prime_mutex);
        }

        mpz_add(candidate, candidate, data->step);
        __sync_fetch_and_add(&primes_checked, 1);

        if (primes_checked % UPDATE_INTERVAL == 0) {
            print_status();
        }
    }

    mpz_clears(candidate, mersenne, NULL);
    return NULL;
}

void print_status(void) {
    time_t current_time = time(NULL);
    double elapsed_time = difftime(current_time, start_time);
    double primes_per_second = primes_checked / elapsed_time;

    printf(ANSI_COLOR_CYAN "\rCurrent n: %llu | " ANSI_COLOR_YELLOW "Primes checked: %llu | " ANSI_COLOR_GREEN "%.2f primes/second" ANSI_COLOR_RESET, 
           current_n, primes_checked, primes_per_second);
    fflush(stdout);
}

int main(int argc, char* argv[]) {
    num_threads = 1;
    unsigned long long initial_n = 3;  // Default starting point (2^3 - 1 = 7 is a Mersenne prime)

    int opt;
    while ((opt = getopt(argc, argv, "t:i:")) != -1) {
        switch (opt) {
            case 't':
                num_threads = atoi(optarg);
                break;
            case 'i':
                initial_n = strtoull(optarg, NULL, 10);
                break;
            default:
                fprintf(stderr, "Usage: %s -t <num_threads> -i <initial_n>\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    if (num_threads < 1 || num_threads > MAX_THREADS) {
        fprintf(stderr, "Number of threads must be between 1 and %d\n", MAX_THREADS);
        exit(EXIT_FAILURE);
    }

    mpz_init(current_prime);
    mpz_set_ui(current_prime, initial_n);
    current_n = initial_n;

    signal(SIGINT, handle_sigint);

    pthread_t threads[MAX_THREADS];
    thread_data_t thread_data[MAX_THREADS];

    start_time = time(NULL);

    for (int i = 0; i < num_threads; i++) {
        mpz_init(thread_data[i].start);
        mpz_init(thread_data[i].step);
        mpz_set_ui(thread_data[i].start, initial_n + i);
        mpz_set_ui(thread_data[i].step, num_threads);
        thread_data[i].thread_id = i;

        if (pthread_create(&threads[i], NULL, find_mersenne_primes, &thread_data[i]) != 0) {
            perror("Failed to create thread");
            exit(EXIT_FAILURE);
        }
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    printf("\n\nSearch completed.\n");

    for (int i = 0; i < num_threads; i++) {
        mpz_clear(thread_data[i].start);
        mpz_clear(thread_data[i].step);
    }

    mpz_clear(current_prime);

    return 0;
}
