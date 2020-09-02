

#ifndef benchmarking_h
#define benchmarking_h

#include <stdio.h>
#include <time.h>
#include "unit_tests.h"

void benchmark_multiplication_bruteforce (void);

void benchmark_multiplication_opt (void);

void benchmark_proof_bruteforce (void);

void benchmark_proof_opt (void);

void benchmark_both_proof (void);



// **********************************************************************
// LIWEI'S METHOD
// **********************************************************************

void benchmark_multiplication_opt_Liwei (void);

void compute_power_2T_wave (mpz_t x, const unsigned long int T, const mpz_t N1, const mpz_t N2, vector* saves, mpz_t out);

void multiplication_wave_trial (void);
void multiplication_wave_trial_while (void);
#endif /* benchmarking_h */
