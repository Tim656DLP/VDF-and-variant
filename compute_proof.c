

#include "compute_proof.h"

//#################################################################################################
void compute_x_i(mpz_t x_in, mpz_t mu_i, mpz_t x_out, mpz_t N, mpz_t r){
    mpz_t res;
    mpz_init(res);
    mpz_powm(res, x_in, r, N);
    mpz_mul(x_out, res, mu_i);
    mpz_mod(x_out, x_out, N);
}

//#################################################################################################
void compute_y_i(mpz_t y_in, mpz_t mu_i, mpz_t y_out, mpz_t N, mpz_t r) {
    mpz_t res;
    mpz_init(res);
    mpz_powm(res, mu_i, r, N);
    mpz_mul(y_out, res, y_in);
    mpz_mod (y_out, y_out, N);
}

//#################################################################################################
void compute_proof_brute_force (const mpz_t x, const mpz_t y, vector* pi, unsigned long int T, mpz_t N) {
    //printf("PROVER \n");
    unsigned long int t = greatest_bit_position(T);
    //printf("T is : %lu \n", T);
    unsigned long int exp = T;
    mpz_t mu;
    mpz_t x_next;
    mpz_t y_next;
    mpz_t exp_mpz;
    mpz_t r;
    mpz_inits(mu, x_next, y_next, r, exp_mpz, NULL);
    mpz_set(x_next, x);
    mpz_set(y_next, y);
    for (int i = 0; i < t; ++ i) {
        mpz_set_ui(r, 0);
        //printf("-------- \n");
        mpz_set_d(exp_mpz, exp);
        exp = exp/2;
        compute_power_2T(x_next, exp, N, mu);
        
        //printf ("mu in proof comp : %lu \n", mpz_get_ui(mu));
        //printf ("exponent in proof comp : %lu \n", mpz_get_ui(exp_mpz));

        hash_function(x_next, exp_mpz, y_next, mu, r);
        compute_x_i(x_next, mu, x_next, N, r);
        compute_y_i(y_next, mu, y_next, N, r);
        //printf("x in ver computation : %lu \n", mpz_get_ui(x_next));
        //printf("y in ver computation : %lu \n", mpz_get_ui(y_next));
        //printf("hash in ver computation : %lu \n", mpz_get_ui(r));

        vector_push(pi, mu);
    }
}

//#################################################################################################
void compute_proof_opt (const mpz_t x, const mpz_t y, vector* pi, unsigned long int T, mpz_t N, vector* saves) {
    
    unsigned long int s = compute_s_parameter(T, 100);
    unsigned long int t = greatest_bit_position(T);

    
    mpz_t mu;
    mpz_t x_next;
    mpz_t y_next;
    mpz_t exp_mu;
    mpz_t exp_x;
    mpz_t r;
    mpz_t T_mpz;
    mpz_t exp_mpz;
    mpz_inits(mu, x_next, y_next, exp_mu, exp_x, r, T_mpz, exp_mpz, NULL);
    mpz_set(x_next, x);
    mpz_set(y_next, y);
    mpz_set_d(exp_x, 1);
    mpz_set_d(exp_mu, 0);
    mpz_set_ui(T_mpz, T);
    
    // computation of the proof with saved results
    unsigned long int next_T = T;
    for (unsigned long int i = 0; i < s; ++i) {
        next_T = next_T/2;
        mpz_setbit(exp_mu, next_T);
        mpz_mul(exp_mu, exp_x, exp_mu);
        exponentiation_for_proof(x, T, N, exp_mu, saves, mu);
        
        //printf ("mu in proof comp : %lu \n", mpz_get_ui(mu));
        //printf ("exponent in proof comp : %lu \n", mpz_get_ui(T_mpz));
        
        hash_function(x_next, T_mpz, y_next, mu, r);
        mpz_mul(exp_x, exp_x, r);
        mpz_add(exp_x, exp_x, exp_mu);
        mpz_set_ui(T_mpz, next_T);
        vector_push(pi, mu);
        compute_x_i(x_next, mu, x_next, N, r);
        compute_y_i(y_next, mu, y_next, N, r);
        
        
        
        //printf("x in ver computation : %lu \n", mpz_get_ui(x_next));
        //printf("y in ver computation : %lu \n", mpz_get_ui(y_next));
        //printf("hash in ver computation : %lu \n", mpz_get_ui(r));
        mpz_set_ui(exp_mu, 0);
        mpz_set_ui(r, 0);
        
    }
    
    // computation of the proof brute force
    for (unsigned long int j = s ; j < t ; ++j) {
        mpz_set_ui(r, 0);
        //printf("-------- \n");
        mpz_set_d(exp_mpz, next_T);
        next_T = next_T/2;
        compute_power_2T(x_next, next_T, N, mu);
        
        //printf ("mu in proof comp : %lu \n", mpz_get_ui(mu));
        //printf ("exponent in proof comp : %lu \n", mpz_get_ui(exp_mpz));
        
        hash_function(x_next, exp_mpz, y_next, mu, r);
        compute_x_i(x_next, mu, x_next, N, r);
        compute_y_i(y_next, mu, y_next, N, r);
        //printf("x in ver computation : %lu \n", mpz_get_ui(x_next));
        //printf("y in ver computation : %lu \n", mpz_get_ui(y_next));
        //printf("hash in ver computation : %lu \n", mpz_get_ui(r));
        
        vector_push(pi, mu);
    }

    
    /*
    //printf("Starting computing proof");
    mpz_t mu;
    mpz_t x_next;
    mpz_t y_next;
    mpz_t exp_mu;
    mpz_t exp_x;
    mpz_t r;
    mpz_t T_mpz;
    mpz_inits(mu, x_next, y_next, exp_mu, exp_x, r, T_mpz, NULL);
    mpz_set(x_next, x);
    mpz_set(y_next, y);
    mpz_set_d(exp_x, 1);
    mpz_set_d(exp_mu, 0);
    mpz_set_ui(T_mpz, T);
    for (unsigned long int i = T/2; i >= 1; i = i/2) {
        //printf("at the %lu th iteration \n", i);
        //printf("-------- \n");
        mpz_setbit(exp_mu, i);
        mpz_mul(exp_mu, exp_x, exp_mu);
        exponentiation_for_proof(x, T, N, exp_mu, saves, mu);
        
        //printf ("mu in proof comp : %lu \n", mpz_get_ui(mu));
        //printf ("exponent in proof comp : %lu \n", mpz_get_ui(T_mpz));
        
        hash_function(x_next, T_mpz, y_next, mu, r);
        mpz_mul(exp_x, exp_x, r);
        mpz_add(exp_x, exp_x, exp_mu);
        mpz_set_ui(T_mpz, i);
        vector_push(pi, mu);
        compute_x_i(x_next, mu, x_next, N, r);
        compute_y_i(y_next, mu, y_next, N, r);
       
        

        //printf("x in ver computation : %lu \n", mpz_get_ui(x_next));
        //printf("y in ver computation : %lu \n", mpz_get_ui(y_next));
        //printf("hash in ver computation : %lu \n", mpz_get_ui(r));
        mpz_set_ui(exp_mu, 0);
        mpz_set_ui(r, 0);
    }
     */
    
    return;
}

//#################################################################################################

