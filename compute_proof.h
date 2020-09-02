

#ifndef compute_proof_h
#define compute_proof_h


#include "hash_function.h"

/*
 * Function:  compute_x_i
 * --------------------
 *  computes the value of x_(i+1) = x^r * mu_i
 *
 *  x_in: the previous occurence of x' (i.e x_(i))
 *  mu_i: the current value of mu_i (i.e mu_(i))
 *  x_out: where the result is to be saved
 *  N: modulo of the class
 *  r: the hash value produced by the hash function
 */
void compute_x_i(mpz_t x_in, mpz_t mu_i, mpz_t x_out, mpz_t N, mpz_t r);

/*
 * Function:  compute_y_i
 * --------------------
 *  computes the value of y_(i+1) = mu_i^r * y_in
 *
 *  y_in: the previous occurence of y' (i.e y_(i))
 *  mu_i: the current value of mu_i (i.e mu_(i))
 *  y_out: where the result is to be saved
 *  N: modulo of the class
 *  r: the hash value produced by the hash function
 */
void compute_y_i(mpz_t y_in, mpz_t mu_i, mpz_t y_out, mpz_t N, mpz_t r);

/*
 * Function:  compute_proof_brute_force
 * --------------------
 *  computes the proof of the VDF without using previously stored result
 *
 *  x_in: the previous occurence of x' (i.e x_(i))
 *  y_in: the previous occurence of y' (i.e y_(i))
 *  pi: vector where the proof is to be saved
 *  T: exponent T/2^i (i.e initialy is equal to T initial)
 *  N: modulo of the class
 */
void compute_proof_brute_force (const mpz_t x, const mpz_t y, vector* pi, unsigned long int T, mpz_t N);

/*
 * Function:  compute_proof_opt
 * --------------------
 *  computes the proof of the VDF without using previously stored result
 *
 *  TO BE COMPLETED
 */
void compute_proof_opt (const mpz_t x, const mpz_t y, vector* pi, unsigned long int T, mpz_t N, vector* saves);


#endif /* compute_proof_h */
