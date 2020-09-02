

#include "helper.h"
#include <math.h>

const unsigned long int STEP = 8;
const unsigned long int lambda = 100;

//#################################################################################################
void compute_exponentiation(mpz_t x, const mpz_t exp, const mpz_t N, mpz_t out) {
    // tests if odd or not ?
    // if T > 1, then recursive operations
    if ( mpz_cmp_ui(exp, 1) > 0 ) {
        mpz_t temp;
        mpz_t new_exp;
        mpz_inits(temp, new_exp, NULL);
        
        // If T is odd, compute using T/2-1
        if(mpz_odd_p(exp) > 0) {
            mpz_sub_ui(new_exp, exp, 1);
            mpz_cdiv_q_ui(new_exp, new_exp, 2);
            compute_exponentiation(x, new_exp, N, temp);
            // temp put allready in modulo
            mpz_mul(out,temp,temp);
            mpz_mul(out, out, x);
            
        } else {
            mpz_cdiv_q_ui(new_exp, exp, 2);
            compute_exponentiation(x, new_exp, N, temp);
            mpz_mul(out, temp,temp);
        }
        // if T = 1, do nothing and simply return y = x
    } else {
        mpz_set(out, x);
    }
}

//#################################################################################################
void compute_power_2T (mpz_t x, const unsigned long int T, const mpz_t N, mpz_t out){
    if (T == 0) {
        mpz_set(out, x);
    } else {
        mpz_t exp;
        mpz_init(exp);
        mpz_setbit(exp, T);
        mpz_powm(out, x, exp, N);
    }
}

//#################################################################################################
unsigned long int compute_s_parameter (unsigned long int T, unsigned long int lambda) {
    int t = (int) log2(T);
    double temp = t/2.0 - log2(t*lambda)/2.0;
    if (temp < 0) {
        return t-3;
    } else {
        return (int) temp;
    }
}

//#################################################################################################
void compute_power_2T_opt (mpz_t x, const unsigned long int T, const mpz_t N, vector* saves, mpz_t out) {
    unsigned long int s = compute_s_parameter(T, lambda);
    unsigned long int power_2STEP = pow(2, s);
    //printf("%d", T/power_2STEP);
    unsigned long int exp_space_between_values = T/power_2STEP;
    
    mpz_t res;
    mpz_init(res);
    compute_power_2T(x, exp_space_between_values, N, res);
    vector_push(saves, res);
    
    for (int i = 0; i < power_2STEP - 1; ++i) {
        compute_power_2T(res, exp_space_between_values, N, res);
        if ( !vector_push(saves, res) ) {
            printf("Unable to save intermediate computation");
        }
    }
    mpz_set(out, res);
}

//#################################################################################################
void exponentiation_for_proof (const mpz_t x, unsigned long int T, mpz_t N, mpz_t exp, vector* saves, mpz_t out) {
    unsigned long int s = compute_s_parameter(T, lambda);
    unsigned long int power_2STEP = pow(2, s);
    unsigned long int exp_space_between_values = T/power_2STEP;
    unsigned long int number_ones = mpz_popcount(exp);
    
    mpz_t exp_copy;
    mpz_init_set(exp_copy, exp);
  
    mpz_t res;
    mpz_init_set_ui(res, 1);
    
    // to be deleted
    //int number_loops = 0;
    
    
    for (unsigned long int i = 0 ; i < number_ones ; ++i) {
       // to be deleted
        //++number_loops;

        unsigned long int q = mpz_sizeinbase(exp_copy, 2) - 1;
        if (q < exp_space_between_values) {
            mpz_t temp;
            mpz_init(temp);
            mpz_powm(temp, x, exp_copy, N);
            mpz_mul(res, res, temp);
            mpz_mod(res, res, N);
            break;
        } else {
            mpz_t saved_value;
            mpz_init(saved_value);
            mpz_clrbit(exp_copy, q);
            //exponentiation_for_proof(x, T, N, exp_copy, saves, rest);
            
            unsigned long int k = q/exp_space_between_values;
            unsigned long int rest_exp = q % exp_space_between_values;
            
            if (k <= power_2STEP) { //|| (k == power_2STEP && mpz_cmp_ui(exp_copy, 0) == 0)) {
                if (k <= saves->size) {
                    vector_get(saves, k - 1, saved_value);
                } else {
                    printf("Index given too big, unable to reach correct saved value. \n");
                }
            } else {
                //printf("computing a greater exponent than pre-computed \n");
                unsigned long int alpha = k/power_2STEP;
                unsigned long int beta = k % power_2STEP;
                
                vector_get(saves, power_2STEP - 1, saved_value);
                // if condition
                compute_power_2T(saved_value, T*(alpha -1), N, saved_value);
                mpz_mod(saved_value, saved_value, N);
                
                // if condition
                compute_power_2T(saved_value, exp_space_between_values * beta, N, saved_value);
                mpz_mod(saved_value, saved_value, N);

                
            }
            // if here as well: choose N1 or N2
            compute_power_2T(saved_value, rest_exp, N, saved_value);
            mpz_mod(saved_value, saved_value, N);

            mpz_mul(res, res, saved_value);
            mpz_mod(res, res, N);
        }

    }
    
    // to be deleted 
    //printf("number of looping : %d \n", number_loops);

    mpz_mod(out, res, N);

    return;

}
/*void exponentiation_for_proof (const mpz_t x, unsigned long int T, mpz_t N, mpz_t exp, vector* saves, mpz_t out) {
    
    unsigned long int power_2STEP = pow(2, STEP);
    unsigned long int exp_space_between_values = T/power_2STEP;
    mpz_t exp_copy;
    mpz_init_set(exp_copy, exp);
    unsigned long int q = mpz_sizeinbase(exp, 2) - 1;

    
    if (q < exp_space_between_values) {
        mpz_powm(out, x, exp, N);
        return;
    } else {
        mpz_t rest;
        mpz_t saved_value;
        mpz_inits(rest, saved_value, NULL);
        mpz_clrbit(exp_copy, q);
        mpz_powm(rest, x, exp_copy, N);
        //exponentiation_for_proof(x, T, N, exp_copy, saves, rest);

        unsigned long int k = q/exp_space_between_values;
        unsigned long int rest_exp = q % exp_space_between_values;

        if (k <= power_2STEP) { //|| (k == power_2STEP && mpz_cmp_ui(exp_copy, 0) == 0)) {
            if (k <= saves->size) { 
                vector_get(saves, k - 1, saved_value);

            } else {
                printf("Index given too big, unable to reach correct saved value. \n");
            }
        } else {
            unsigned long int alpha = k/power_2STEP;
            unsigned long int beta = k % power_2STEP;

            vector_get(saves, power_2STEP - 1, saved_value);
            compute_power_2T(saved_value, T*(alpha -1), N, saved_value);
            compute_power_2T(saved_value, exp_space_between_values * beta, N, saved_value);

        }
        compute_power_2T(saved_value, rest_exp, N, saved_value);

        mpz_mul(out, rest, saved_value);
        mpz_mod(out, out, N);
        return;
    }
    
}*/

//#################################################################################################
int check_quatratic_residue(const mpz_t var, const mpz_t N) {
    if( (var < 0) && (mpz_jacobi(var, N) != 1) )  {
        return 1;

    } else {
       return 0;
    }
}

//#################################################################################################
unsigned long int greatest_bit_position(unsigned long int num) {
    if (num == 0) {
        printf("Num is zero, greatest bit doesn't exists");
        return 0;
    }
    unsigned long int m;
    m = num;
    m = m | m >> 1;
    m = m | m >> 2;
    m = m | m >> 4;
    m = m | m >> 8;
    m = m | m >> 16;
    m = m | m >> 32;
    m = m & ((~m >> 1)); //^0x8000000000000000);

    unsigned long int res = 0;
    
    while ( m != 0) {
        m = m >> 1;
        res++;
    }
    return res - 1 ;

}

void generate_RSAmodulus (unsigned long int size, mpz_t out) {
    mpz_t safeprime_1;
    mpz_t safeprime_2;
    mpz_t random_number;
    mpz_inits(random_number, safeprime_1, safeprime_2, NULL);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    int count = 0;
    
    while (1) {
        mpz_urandomb(random_number, state, size);
        //printf("random number under analysis is : %lu \n", mpz_get_ui(random_number));
        
        int prob_prime = mpz_probab_prime_p(random_number, 20);
        if (prob_prime > 0) {
            mpz_t potential_safeprime;
            mpz_init(potential_safeprime);
            mpz_mul_ui(potential_safeprime, random_number, 2);
            mpz_add_ui(potential_safeprime, potential_safeprime, 1);
            int prob_safe = mpz_probab_prime_p(potential_safeprime, 20);
            if (prob_safe > 0) {
                if (count == 0) {
                    mpz_set(safeprime_1, potential_safeprime);
                    //gmp_printf ("Here is a safe prime is an mpz %Zd\n",safeprime_1);
                    count += 1;
                } else {
                    mpz_set(safeprime_2, potential_safeprime);
                    //gmp_printf ("Here is a safe prime is an mpz %Zd\n",safeprime_2);
                    mpz_mul(out, safeprime_1, safeprime_2);
                    gmp_printf ("Here is the prime modulus %Zd\n",out);
                    return;
                }
            }
        }
    }
}

void generate_safeprime(unsigned long int size, mpz_t out){
   
    mpz_t random_number;
    mpz_init(random_number);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    
    while (1) {
        mpz_urandomb(random_number, state, size);
        printf("random number under analysis is : %lu \n", mpz_get_ui(random_number));
        
        int prob_prime = mpz_probab_prime_p(random_number, 20);
        if (prob_prime > 0) {
            mpz_t potential_safeprime;
            mpz_init(potential_safeprime);
            mpz_mul_ui(potential_safeprime, random_number, 2);
            mpz_add_ui(potential_safeprime, potential_safeprime, 1);
            int prob_safe = mpz_probab_prime_p(potential_safeprime, 20);
            if (prob_safe > 0) {
                mpz_set(out, potential_safeprime);
                gmp_printf ("Here is a safe prime is an mpz %Zd\n",out);
                return;
            }
        }
    }
    
    //##################################################
    /*
    mpz_t greater;
    mpz_inits(greater, NULL);
    mpz_setbit(greater, size);
    
    while (1) {
        //printf("Try to find a prime \n");
        mpz_nextprime(out, greater);
        
        mpz_t power;
        mpz_t exp;
        mpz_t mod;
        mpz_inits(power, exp, mod, NULL);
        
        mpz_mul_ui(exp, out, 2);
        mpz_set_ui(power, 2);
        mpz_mul_ui(mod, out, 2);
        mpz_add_ui(mod, mod, 1);
        mpz_powm(power, power, exp, mod);
    
        if (mpz_cmp_ui(power, 1) == 0) {
            mpz_set(out, mod);
            break;
        }
        mpz_set(greater, out);
    }
    
    gmp_printf ("%s is an mpz %Zd\n", "here", out);
     */
    
    //##################################################
    /*
    mpz_t greater;
    mpz_inits(greater, NULL);
    mpz_setbit(greater, size);
    
    while (1) {
        //printf("Try to find a prime \n");
        mpz_nextprime(out, greater);
        
        //printf("Prime under study : %lu \n", mpz_get_ui(out));
        mpz_t potential_safeprime;
        mpz_init(potential_safeprime);

        mpz_sub_ui(potential_safeprime, out, 1);
        mpz_div_ui(potential_safeprime, potential_safeprime, 2);
        //printf("Potential safeprime : %lu", mpz_get_ui(potential_safeprime));
        int res = 16;
        int prob = mpz_probab_prime_p(potential_safeprime, res);
        
        //printf("Prob. is : %d \n", prob);
        if (prob >= 1) {
            break;
        }
        mpz_set(greater, out);
        //printf("Did not found a prime. Trying again \n");
    }
    
    gmp_printf ("%s is an mpz %Zd\n", "here", out);
    */
    
    //##################################################
    /*
    mpz_t greater;
    mpz_inits(greater, NULL);
    mpz_setbit(greater, size);
    
    while (1) {
        printf("Try to find a prime \n");
        mpz_nextprime(out, greater);
        mpz_t potential_safeprime;
        mpz_init(potential_safeprime);
        mpz_mul_ui(potential_safeprime, out, 2);
        mpz_sub_ui(potential_safeprime, potential_safeprime, 1);
        int prob = mpz_probab_prime_p(potential_safeprime, 16);
        if (prob > 0) {
            mpz_set(out, potential_safeprime);
            break;
        }
        printf("Did not found a prime. Trying again \n");
    }
    
    printf("safe prime found %lu : ", mpz_get_ui(out));
    */
    
    return;
}

// ***************************************************************************
// LIWEI'S METHOD
// ***************************************************************************


//#################################################################################################
void compute_power_2T_opt_Liwei (mpz_t x, const unsigned long int T, const mpz_t N1, const mpz_t N2, vector* saves, mpz_t out) {
    unsigned long int s = compute_s_parameter(T, lambda);
    unsigned long int power_2STEP = pow(2, s);
    //printf("%d", T/power_2STEP);
    unsigned long int exp_space_between_values = T/power_2STEP;
    
    mpz_t res;
    mpz_init(res);
    compute_power_2T(x, exp_space_between_values, N1, res);
    vector_push(saves, res);
    
    int counter = 1;
    for (int i = 0; i < power_2STEP - 1; ++i) {
        if (counter == 1){
            compute_power_2T(res, exp_space_between_values, N2, res);
            counter = 0;
        } else {
            compute_power_2T(res, exp_space_between_values, N1, res);
            counter = 1;
        }
        if ( !vector_push(saves, res) ) {
            printf("Unable to save intermediate computation");
        }
    }
    mpz_set(out, res);
}



