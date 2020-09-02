

#ifndef helper_h
#define helper_h

#include <string.h>
#include <math.h>

#include "vector.h"


/*
 * Function:  compute_exponentiation
 * --------------------
 * computes recusively computes x^exp
 *
 *  x: the variable to be raised to the exponent
 *  exp: the exponent
 *  N: the modulo of the space in which computation are done
 *  out: the variable in which the result is saved
 */
void compute_exponentiation(mpz_t x, const mpz_t exp, const mpz_t N, mpz_t out);


/*
 * Function:  compute_power_2T
 * --------------------
 * computes recusively out = x^(2^T)
 *
 *  x: the number to be raised to the power of 2^T
 *  T: the sub-exponential power
 *  N: the modulo of the space in which computation are done
 *  out: the variable in which the result is saved
 */
void compute_power_2T(mpz_t x, const unsigned long int T, const mpz_t N, mpz_t out);

/*
 * Function:  compute_s_parameter
 * --------------------
 * compute the s parameter of the proof. As a reminder pow(2, s) is the number of saved values
 * assumption : T is less than or equal to 64 bits & size_t is encoded on 64 bits
 *
 *  T: the variable used to compute x^(2^T) in the Time lock puzzle
 *  lambda : the security parameter
 *  exp: the exponent
 *  saves: the vector in which the intermediate result are saved (note : first value of the vector is for exponent T/(2^step), then 2*T/(2^step),
 *
 *  returns the s parameter as an unsigned long int
 */
unsigned long int compute_s_parameter (unsigned long int T, unsigned long int lambda);

/*
 * Function:  compute_power_2T_opt
 * --------------------
 * computes recusively out = x^(2^T) and saves some results into the vector saves
 *
 *  x: the number to be raised to the power of 2^T
 *  T: the sub-exponential power
 *  N: the modulo of the space in which computation are done
 *  saves: the vector in which the intermediate result are saved (note : first value of the vector is for exponent T/(2^step), then 2*T/(2^step),
 *  out: the variable in which the result is saved
 */
void compute_power_2T_opt(mpz_t x, const unsigned long int T, const mpz_t N, vector* saves, mpz_t out);


/*
 * Function:  exponentiation_for_proof
 * --------------------
 * computes x^(exp) using previously saved result
 * assumption : T is less than or equal to 64 bits & size_t is encoded on 64 bits
 *
 *  x: the number to be raised to the power of exp
 *  T: the variable used to compute x^(2^T) in the Time lock puzzle
 *  N: the modulo of the space in which computation are done
 *  exp: the exponent
 *  saves: the vector in which the intermediate result are saved (note : first value of the vector is for exponent T/(2^step), then 2*T/(2^step),
 *
 */
void exponentiation_for_proof (const mpz_t x, unsigned long int T, mpz_t N, mpz_t exp, vector* saves, mpz_t out);

/*
 * Function:  check_quatratic_residue
 * --------------------
 * checks that the given variable var is in the positive QR_N (done by computing the Jacobi symbols of the variable)
 *
 *  var: the number to be checked if in positve QR_N
 *  N: the modulo of the group
 *
 *  returns: 0 if var is in positive QR_N else 1
 */
int check_quatratic_residue(const mpz_t var, const mpz_t N);

/*
 * Function:  greatest_bit_position
 * --------------------
 * computes the position of the most significant bit of num
 *
 *  num : an unsigned long int variable (encoded on 64 bits)
 *
 *  returns: the positon of the most significant bit 
 */
unsigned long int greatest_bit_position(unsigned long int num);

/*
 * Function:  generate_RSAmodulus
 * --------------------
 * computes the RSA modulus N of the VDF (i.e find two safeprime and return
 * their product
 *
 *  size : size of each safe prime to be found. The final size of the modulus is size*2
 *
 *  returns: an RSA modulus of size size*2
 */
void generate_RSAmodulus (unsigned long int size, mpz_t out);

void generate_safeprime(unsigned long int size, mpz_t out);

// **********************************************************************
// LIWEI'S METHOD
// **********************************************************************

void compute_power_2T_opt_Liwei (mpz_t x, const unsigned long int T, const mpz_t N1, const mpz_t N2, vector* saves, mpz_t out);

#endif /* helper_h */
