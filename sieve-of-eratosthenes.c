#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>
#include <string.h>

#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_RESET "\x1b[0m"

#define NUM_THREADS 4

typedef struct {
    unsigned long long start;
    unsigned long long end;
    bool *sieve;
    unsigned long long limit;
} ThreadData;

unsigned long long primes_found = 0;
pthread_mutex_t primes_mutex = PTHREAD_MUTEX_INITIALIZER;

void *sieve_thread(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    unsigned long long i, j;

    for (i = data->start; i <= data->end && i * i <= data->limit; i++) {
        if (!data->sieve[i]) {
            for (j = i * i; j <= data->limit; j += i) {
                data->sieve[j] = true;
            }
        }
    }

    return NULL;
}

void count_primes(bool *sieve, unsigned long long limit) {
    unsigned long long i;
    unsigned long long count = 0;

    for (i = 2; i <= limit; i++) {
        if (!sieve[i]) {
            count++;
        }
    }

    pthread_mutex_lock(&primes_mutex);
    primes_found = count;
    pthread_mutex_unlock(&primes_mutex);
}

int main() {
    unsigned long long limit;
    printf("Enter the upper limit for prime number search: ");
    scanf("%llu", &limit);

    bool *sieve = calloc(limit + 1, sizeof(bool));
    if (!sieve) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    sieve[0] = sieve[1] = true;

    pthread_t threads[NUM_THREADS];
    ThreadData thread_data[NUM_THREADS];

    unsigned long long chunk_size = (limit - 1) / NUM_THREADS;
    clock_t start_time = clock();
    time_t last_update = time(NULL);

    for (int i = 0; i < NUM_THREADS; i++) {
        thread_data[i].start = 2 + i * chunk_size;
        thread_data[i].end = (i == NUM_THREADS - 1) ? limit : (i + 1) * chunk_size + 1;
        thread_data[i].sieve = sieve;
        thread_data[i].limit = limit;

        if (pthread_create(&threads[i], NULL, sieve_thread, &thread_data[i]) != 0) {
            fprintf(stderr, "Failed to create thread %d\n", i);
            free(sieve);
            return 1;
        }
    }

    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    count_primes(sieve, limit);

    clock_t end_time = clock();
    double cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

    printf("\n" ANSI_COLOR_CYAN "Prime numbers up to %llu:" ANSI_COLOR_RESET "\n", limit);
    unsigned long long printed = 0;
    for (unsigned long long i = 2; i <= limit; i++) {
        if (!sieve[i]) {
            printf("%llu ", i);
            printed++;
            if (printed % 10 == 0) printf("\n");
        }
    }

    printf("\n\n" ANSI_COLOR_YELLOW "Total prime numbers found: %llu" ANSI_COLOR_RESET "\n", primes_found);
    printf(ANSI_COLOR_GREEN "Time taken: %.2f seconds (mostly acurate most of the time except when it is not)" ANSI_COLOR_RESET "\n", cpu_time_used);
    printf(ANSI_COLOR_CYAN "Primes per second: %.2f" ANSI_COLOR_RESET "\n", primes_found / cpu_time_used);

    free(sieve);
    return 0;
}
