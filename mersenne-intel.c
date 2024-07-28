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
#include <sched.h>
#include <numa.h>
#include <atomic>
#include <x86intrin.h>

#define TOP_LEVEL_THREADS 16
#define WORKER_THREADS_PER_TOP 24
#define TOTAL_THREADS (TOP_LEVEL_THREADS * WORKER_THREADS_PER_TOP)
#define MILLER_RABIN_ITERATIONS 40
#define UPDATE_INTERVAL 10000
#define CHUNK_SIZE 1000000

// ANSI color codes (unchanged)

typedef struct {
    mpz_t start;
    mpz_t step;
    int thread_id;
    int numa_node;
    pthread_t* worker_threads;
    std::atomic<unsigned long long>* local_primes_checked;
} thread_data_t;

volatile sig_atomic_t keep_running = 1;
pthread_mutex_t prime_mutex = PTHREAD_MUTEX_INITIALIZER;
mpz_t current_prime;
std::atomic<unsigned long long> current_n;
int num_top_threads;
std::atomic<unsigned long long> primes_checked{0};
time_t start_time;

// Function prototypes
void print_status(void);
void handle_sigint(int sig);
int miller_rabin(mpz_t n, int iterations);
void* find_mersenne_primes_worker(void* arg);
void* find_mersenne_primes_controller(void* arg);

// Optimized Miller-Rabin primality test
int miller_rabin(mpz_t n, int iterations) {
    if (mpz_cmp_ui(n, 2) < 0) return 0;
    if (mpz_cmp_ui(n, 2) == 0) return 1;
    if (mpz_even_p(n)) return 0;

    mpz_t a, y, n_minus_one;
    mpz_inits(a, y, n_minus_one, NULL);

    unsigned long int s = 0;
    mpz_sub_ui(n_minus_one, n, 1);
    mpz_t r;
    mpz_init_set(r, n_minus_one);

    while (mpz_even_p(r)) {
        mpz_fdiv_q_2exp(r, r, 1);
        s++;
    }

    gmp_randstate_t rnd;
    gmp_randinit_mt(rnd);
    gmp_randseed_ui(rnd, _rdtsc());

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

    mpz_clears(a, y, r, n_minus_one, NULL);
    gmp_randclear(rnd);
    return is_prime;
}

void* find_mersenne_primes_worker(void* arg) {
    thread_data_t* data = (thread_data_t*)arg;
    mpz_t candidate, mersenne;
    mpz_inits(candidate, mersenne, NULL);
    mpz_set(candidate, data->start);

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(data->thread_id, &cpuset);
    pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);

    unsigned long long local_checked = 0;

    while (keep_running) {
        for (int i = 0; i < CHUNK_SIZE && keep_running; i++) {
            mpz_ui_pow_ui(mersenne, 2, mpz_get_ui(candidate));
            mpz_sub_ui(mersenne, mersenne, 1);

            if (miller_rabin(mersenne, MILLER_RABIN_ITERATIONS)) {
                pthread_mutex_lock(&prime_mutex);
                if (mpz_cmp(candidate, current_prime) > 0) {
                    mpz_set(current_prime, candidate);
                    current_n.store(mpz_get_ui(candidate));
                    printf(ANSI_COLOR_GREEN "\nFound Mersenne prime: 2^%llu - 1\n" ANSI_COLOR_RESET, current_n.load());
                    mpz_out_str(stdout, 10, mersenne);
                    printf("\n");
                }
                pthread_mutex_unlock(&prime_mutex);
            }

            mpz_add(candidate, candidate, data->step);
            local_checked++;
        }

        data->local_primes_checked->fetch_add(local_checked, std::memory_order_relaxed);
        local_checked = 0;

        if (data->local_primes_checked->load(std::memory_order_relaxed) % UPDATE_INTERVAL == 0) {
            primes_checked.fetch_add(data->local_primes_checked->exchange(0, std::memory_order_relaxed), std::memory_order_relaxed);
            print_status();
        }
    }

    mpz_clears(candidate, mersenne, NULL);
    return NULL;
}

void* find_mersenne_primes_controller(void* arg) {
    thread_data_t* data = (thread_data_t*)arg;
    
    // Set NUMA affinity for the controller thread
    numa_run_on_node(data->numa_node);

    // Create worker threads
    for (int i = 0; i < WORKER_THREADS_PER_TOP; i++) {
        thread_data_t* worker_data = (thread_data_t*)malloc(sizeof(thread_data_t));
        memcpy(worker_data, data, sizeof(thread_data_t));
        worker_data->thread_id = data->thread_id * WORKER_THREADS_PER_TOP + i;
        mpz_add_ui(worker_data->start, worker_data->start, i * TOP_LEVEL_THREADS);
        mpz_mul_ui(worker_data->step, worker_data->step, WORKER_THREADS_PER_TOP);

        if (pthread_create(&data->worker_threads[i], NULL, find_mersenne_primes_worker, worker_data) != 0) {
            perror("Failed to create worker thread");
            exit(EXIT_FAILURE);
        }
    }

    // Wait for worker threads to complete
    for (int i = 0; i < WORKER_THREADS_PER_TOP; i++) {
        pthread_join(data->worker_threads[i], NULL);
    }

    return NULL;
}

// Main function (with slight modifications)
int main(int argc, char* argv[]) {
    num_top_threads = TOP_LEVEL_THREADS;
    unsigned long long initial_n = 3;

    // Parse command-line arguments (unchanged)

    mpz_init(current_prime);
    mpz_set_ui(current_prime, initial_n);
    current_n.store(initial_n);

    signal(SIGINT, handle_sigint);

    pthread_t top_threads[TOP_LEVEL_THREADS];
    thread_data_t thread_data[TOP_LEVEL_THREADS];

    start_time = time(NULL);

    for (int i = 0; i < num_top_threads; i++) {
        mpz_init(thread_data[i].start);
        mpz_init(thread_data[i].step);
        mpz_set_ui(thread_data[i].start, initial_n + i);
        mpz_set_ui(thread_data[i].step, num_top_threads);
        thread_data[i].thread_id = i;
        thread_data[i].numa_node = i % numa_num_configured_nodes();
        thread_data[i].worker_threads = (pthread_t*)malloc(WORKER_THREADS_PER_TOP * sizeof(pthread_t));
        thread_data[i].local_primes_checked = new std::atomic<unsigned long long>(0);

        if (pthread_create(&top_threads[i], NULL, find_mersenne_primes_controller, &thread_data[i]) != 0) {
            perror("Failed to create top-level thread");
            exit(EXIT_FAILURE);
        }
    }

    for (int i = 0; i < num_top_threads; i++) {
        pthread_join(top_threads[i], NULL);
    }

    printf("\n\nSearch completed.\n");

    // Clean up resources
    for (int i = 0; i < num_top_threads; i++) {
        mpz_clear(thread_data[i].start);
        mpz_clear(thread_data[i].step);
        free(thread_data[i].worker_threads);
        delete thread_data[i].local_primes_checked;
    }

    mpz_clear(current_prime);

    return 0;
}

// Print status function (unchanged)
void print_status(void) {
    time_t current_time = time(NULL);
    double elapsed_time = difftime(current_time, start_time);
    double primes_per_second = primes_checked.load() / elapsed_time;

    printf(ANSI_COLOR_CYAN "\rCurrent n: %llu | " ANSI_COLOR_YELLOW "Primes checked: %llu | " ANSI_COLOR_GREEN "%.2f primes/second" ANSI_COLOR_RESET, 
           current_n.load(), primes_checked.load(), primes_per_second);
    fflush(stdout);
}

// Signal handler (unchanged)
void handle_sigint(int sig) {
    keep_running = 0;
}
