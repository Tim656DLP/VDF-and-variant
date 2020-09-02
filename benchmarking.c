

#include "benchmarking.h"

void benchmark_multiplication_bruteforce(void) {
    mpz_t x;
    mpz_t N;
    mpz_t out;
    mpz_inits(x, N, out, NULL);
    mpz_set_ui(x, 450003);
    mpz_set_str(N, "622508016097411244439761737392416325553332376060386097641057707287735055303924775370444497742993919891589684328516458697581692371705394041823221085279874642396220161984177201792600246913045294001429059109375816432086923550995143584650550477", 10);
    
    // creation of the csv file to save result
    char * filename = "bruteforce_multiplication.csv";
    printf("Creating %s file \n",filename);
    FILE *fp;
    
    fp=fopen(filename,"w+");
    clock_t time;

    fprintf(fp,"sub_exponent,time");
    for (int t = 20; t < 21; t++){
        unsigned long int T = pow(2, t);
        fprintf(fp,"\n%lu",T);
        
        time = clock();
        compute_power_2T(x, T, N, out);
        time = clock() - time;
        printf("time is %f\n", (double) time/CLOCKS_PER_SEC);

        fprintf(fp,",%lu",time);
        
    }
    fclose(fp);
    
    printf("%s file created \n",filename);
}

void benchmark_multiplication_opt(void) {
    mpz_t x;
    mpz_t N;
    mpz_t out;
    mpz_inits(x, N, out, NULL);
    mpz_set_ui(x, 350003);
    
    mpz_set_str(N, "622508016097411244439761737392416325553332376060386097641057707287735055303924775370444497742993919891589684328516458697581692371705394041823221085279874642396220161984177201792600246913045294001429059109375816432086923550995143584650550477", 10);
    
    // creation of the csv file to save result
    char * filename = "opt_multiplication.csv";
    printf("Creating %s file \n",filename);
    FILE *fp;
    
    fp=fopen(filename,"w+");
    
    fprintf(fp,"sub_exponent, time");
    for (int t = 13; t < 30; t++){
        unsigned long int T = pow(2, t);
        fprintf(fp,"\n%lu",T);
        vector saves;
        construct_vector(&saves);
        
        clock_t time;
        time = clock();
        compute_power_2T_opt(x, T, N, &saves, out);
        time = clock() - time;
        printf("time is %lu \n", time);

        double time_taken = ((double)time);
        
        fprintf(fp,",%f ",time_taken);
        
    }
    fclose(fp);

    
    printf("%s file created \n",filename);
}

void benchmark_proof_bruteforce (void) {
    mpz_t x;
    mpz_t N;
    mpz_t y;
    mpz_inits(x, N, y, NULL);
    mpz_set_ui(x, 123456);
    mpz_set_str(N, "16158503035655503650357438344334975980222051334857742016065172713762327569433945446598600705761456731844358980460949009747059779575245460547544076193224141560315438683650498045875098875194826053398028819192033784138396109321309878080919047169238085235290822926018152521443787945770532904303776199561965642282844350721009131075592574364243695880618860889291860118718196676371399528079340206559035126612674866232053943452084332329067399812414572508105263197679145256552814549607319534463349117382246918962235483023376108551954504056221074590325455782931366628146267643957799255048569835849362555689794934167531092723509", 10);
    
    // creation of the csv file to save result
    char * filename_opt = "opt_proof.csv";
    printf("Creating %s file \n",filename_opt);
    FILE *fp_opt;
    
    fp_opt=fopen(filename_opt,"w+");
    
    fprintf(fp_opt,"T,time");
    
    // creation of the csv file to save time to compute challenge
    char * filename_mul = "multiplication.csv";
    printf("Creating %s file \n",filename_mul);
    FILE *fp_mul;
    
    fp_mul=fopen(filename_mul,"w+");
    fprintf(fp_mul,"T,time");
    
    // creation of the csv file to compute the verifiction
    char * filename_ver = "verification.csv";
    printf("Creating %s file \n",filename_ver);
    FILE *fp_ver;
    
    fp_ver=fopen(filename_ver,"w+");
    fprintf(fp_ver,"T,time");
    
    
    for (int t = 13; t < 30; t++){
        printf("at iteration : %d \n", t);
        unsigned long int T = pow(2, t);
        vector saves;
        construct_vector(&saves);
        
        // Compute mul
        fprintf(fp_mul,"\n%lu",T);
        
        clock_t time_mul;
        time_mul = clock();
        compute_power_2T_opt(x, T, N, &saves, y);
        time_mul = clock() - time_mul;
        
        fprintf(fp_mul,",%f ", (double) time_mul/CLOCKS_PER_SEC);
        
        
        
        mpz_mod (y, y, N);
        vector pi_opt;
        construct_vector(&pi_opt);
        
        // Compute proof optimized version
        fprintf(fp_opt,"\n%lu",T);
        
        clock_t time_opt;
        time_opt = clock();
        compute_proof_opt(x, y, &pi_opt, T, N, &saves);
        time_opt = clock() - time_opt;
        
        fprintf(fp_opt,",%f ", (double) time_opt/CLOCKS_PER_SEC);
        
        // verify the proof
        clock_t time_ver;
        time_ver = clock();
        compute_proof_opt(x, y, &pi_opt, T, N, &saves);
        time_ver = clock() - time_ver;
        
        fprintf(fp_opt,",%f ", (double) time_opt/CLOCKS_PER_SEC);
        
        
        
        
        
    }
    fclose(fp_opt);
    fclose(fp_mul);
    
    printf("%s file created \n",filename_opt);
    printf("%s file created \n",filename_mul);

}

void benchmark_both_proof (void) {
    mpz_t x;
    mpz_t N;
    mpz_t y;
    mpz_inits(x, N, y, NULL);
    mpz_set_ui(x, 123456);
    mpz_set_str(N, "16158503035655503650357438344334975980222051334857742016065172713762327569433945446598600705761456731844358980460949009747059779575245460547544076193224141560315438683650498045875098875194826053398028819192033784138396109321309878080919047169238085235290822926018152521443787945770532904303776199561965642282844350721009131075592574364243695880618860889291860118718196676371399528079340206559035126612674866232053943452084332329067399812414572508105263197679145256552814549607319534463349117382246918962235483023376108551954504056221074590325455782931366628146267643957799255048569835849362555689794934167531092723509", 10);
    
    // creation of the csv file to save result
    char * filename_opt = "opt_proof.csv";
    printf("Creating %s file \n",filename_opt);
    FILE *fp_opt;
    
    fp_opt=fopen(filename_opt,"w+");
    
    fprintf(fp_opt,"T,time");
    
    char * filename_brute = "bruteforce_proof.csv";
    printf("Creating %s file \n",filename_brute);
    FILE *fp_brute;
    
    fp_brute=fopen(filename_brute,"w+");
    
    fprintf(fp_brute,"T,time");
    
    char * filename_mul = "multiplication.csv";
    printf("Creating %s file \n",filename_mul);
    FILE *fp_mul;
    
    fp_mul=fopen(filename_mul,"w+");
    fprintf(fp_mul,"T,time");

    
    
    for (int t = 13; t < 30; t++){
        printf("at iteration : %d \n", t);
        unsigned long int T = pow(2, t);
        vector saves;
        construct_vector(&saves);
        
        // Compute mul
        fprintf(fp_mul,"\n%lu",T);
        
        clock_t time_mul;
        time_mul = clock();
        compute_power_2T_opt(x, T, N, &saves, y);
        time_mul = clock() - time_mul;
        
        fprintf(fp_mul,",%f ", (double) time_mul/CLOCKS_PER_SEC);

        
        
        mpz_mod (y, y, N);
        vector pi_opt;
        construct_vector(&pi_opt);
        vector pi_brute;
        construct_vector(&pi_brute);
        
        // Compute proof optimized version
        fprintf(fp_opt,"\n%lu",T);
        
        clock_t time_opt;
        time_opt = clock();
        compute_proof_opt(x, y, &pi_opt, T, N, &saves);
        time_opt = clock() - time_opt;
        
        fprintf(fp_opt,",%f ", (double) time_opt/CLOCKS_PER_SEC);
        
        // Compute proof brute force version
        fprintf(fp_brute,"\n%lu",T);

        clock_t time_brute;
        time_brute = clock();
        compute_proof_brute_force(x, y, &pi_brute, T, N);
        time_brute = clock() - time_brute;
        
        fprintf(fp_brute,",%f ", (double) time_brute/CLOCKS_PER_SEC);

        
    }
    fclose(fp_opt);
    fclose(fp_brute);
    
    printf("%s file created \n",filename_opt);
    printf("%s file created \n",filename_brute);

}


// **********************************************************************
// LIWEI'S METHOD
// **********************************************************************

/*
void benchmark_multiplication_opt_Liwei (void) {
    mpz_t x;
    mpz_t N1;
    mpz_t N2;
    mpz_t out;
    mpz_inits(x, N1, N2, out, NULL);
    mpz_set_ui(x, 350003);
    
    mpz_set_str(N1, "14519772012476423465891982448492203124161614999261905964890938397438520474358133001222396013468252639273497391014076774465353149747064046378663658364877682332792216936398941896675303536904501882196414133904147973257160025989192740039367903591217360707468865930770369608248883163613941747293380289534338007988556642276930460260647702107334829120900171946977524054284673452269133833025974184893751704427960235064862223990436511490270330371095192617976695419601386221160501773650442727691167641401336546436476298872397765696993333642940661853695364626096606074201498383021755685061873527964530976679135850914734022134269", 10);
    mpz_set_str(N2, "16586590687361643587612739050606246075883818593688125833577699316220127977080168588304354324356445963624814577300868381401214405858367204503363839369063982428691566418281121188688754598644319796441837767922284499546224705478663209861083006852760229926776853038094754552070328034104193400893245372549030298792777754250648573691618715924791001499873681565005082286508541090268441620053266266867900385747199730993494738607288416717003740955088644237161192943830787190330322349014709473604906447330175421589118526590859764383851683929006141834660086562729323700925281924198630248699094612603939433824336106029747514765083", 10);
    
    // creation of the csv file to save result
    char * filename = "opt_multiplication_LIWEI.csv";
    printf("Creating %s file \n",filename);
    FILE *fp;
    
    fp=fopen(filename,"w+");
    
    fprintf(fp,"sub_exponent, time");
    for (int t = 13; t < 30; t++){
        unsigned long int T = pow(2, t);
        fprintf(fp,"\n%lu",T);
        vector saves;
        construct_vector(&saves);
        
        clock_t time;
        time = clock();
        compute_power_2T_opt_Liwei(x, T, N1, N2, &saves, out);
        time = clock() - time;
        printf("time is %lu \n", time);
        
        double time_taken = ((double)time);
        
        fprintf(fp,",%f ",time_taken);
        
    }
    fclose(fp);
    
    
    printf("%s file created \n",filename);
}
 */


void compute_power_2T_wave (mpz_t x, const unsigned long int T, const mpz_t N1, const mpz_t N2, vector* saves, mpz_t out) {
    
    char * filename = "opt_multiplication_WAVE.csv";
    printf("Creating %s file \n",filename);
    FILE *fp;
    
    fp=fopen(filename,"w+");
    fprintf(fp,"exp, res");

    
    mpz_t res;
    mpz_init(res);
    mpz_set(res, x);
    
    int counter = 0;
    for (int i = 1; i <= T; i++){
        if (counter == 1){
            compute_power_2T(res, 2, N2, res);
            counter = 0;
        } else {
            // res = x
            // T=1: res = res^2 mod N1
            // T=2: res = res^2 mod N2 ...
        
            compute_power_2T(res, 2, N1, res);
            counter = 1;
        }
        char * s = mpz_get_str(NULL, 10, res);
        fprintf(fp, "\n%d, %s", i, s);
        //vector_get(vector* vect,size_t pos, mpz_t out);
        mpz_t temp;
        mpz_init(temp);

        int break_ = 0;
        for (int j = 0; j < saves->allocated; j++) {
            vector_get(saves, j, temp);
            if (!mpz_cmp(temp, res)) {
                printf("loop: %d", i);
                break_ = 1;
                break;
            }
        }
        
        if (break_) {
            break;
        }
        if ( !vector_push(saves, res) ) {
            printf("Unable to save intermediate computation");
        }
        
    }
    mpz_set(out, res);
    fclose(fp);
    

    printf("%s file created \n",filename);
}

void compute_power_2T_wave_while (mpz_t x, const mpz_t N1, const mpz_t N2, vector* saves, mpz_t out) {
    
    char * filename = "opt_multiplication_WAVE_while.csv";
    printf("Creating %s file \n",filename);
    FILE *fp;
    
    fp=fopen(filename,"w+");
    fprintf(fp,"exp, res");
    
    
    mpz_t res;
    mpz_init(res);
    mpz_set(res, x);
    
    int counter = 0;
    int break_ = 0;
    int i = 1;
    while (!break_){
        printf("starting: %d \n", i);
        if (counter == 1){
            compute_power_2T(res, 2, N2, res);
            counter = 0;
        } else {
            // res = x
            // T=1: res = res^2 mod N1
            // T=2: res = res^2 mod N2 ...
            
            compute_power_2T(res, 2, N1, res);
            counter = 1;
        }
        char * s = mpz_get_str(NULL, 10, res);
        fprintf(fp, "\n%d, %s", i, s);
        //vector_get(vector* vect,size_t pos, mpz_t out);
        mpz_t temp;
        mpz_init(temp);
        
        for (int j = 0; j < saves->allocated; j++) {
            vector_get(saves, j, temp);
            if (!mpz_cmp(temp, res)) {
                printf("loop: %d", i);
                break_ = 1;
                break;
            }
        }
        
        if (break_) {
            break;
        }
        if ( !vector_push(saves, res) ) {
            printf("Unable to save intermediate computation");
        }
        i += 1;
    }
    mpz_set(out, res);
    fclose(fp);
    
    
    printf("%s file created \n",filename);
}


void multiplication_wave_trial (void){
    mpz_t x;
    mpz_t N1;
    mpz_t N2;
    mpz_t out;
    mpz_inits(x, N1, N2, out, NULL);
    vector saves;
    construct_vector(&saves);
    
    // CHANGE PARAMETERS HERE
    mpz_set_ui(x, 7);
    unsigned long int T = pow(2, 14);
    mpz_set_str(N1, "14519772012476423465891982448492203124161614999261905964890938397438520474358133001222396013468252639273497391014076774465353149747064046378663658364877682332792216936398941896675303536904501882196414133904147973257160025989192740039367903591217360707468865930770369608248883163613941747293380289534338007988556642276930460260647702107334829120900171946977524054284673452269133833025974184893751704427960235064862223990436511490270330371095192617976695419601386221160501773650442727691167641401336546436476298872397765696993333642940661853695364626096606074201498383021755685061873527964530976679135850914734022134269", 10);
    mpz_set_str(N2, "16586590687361643587612739050606246075883818593688125833577699316220127977080168588304354324356445963624814577300868381401214405858367204503363839369063982428691566418281121188688754598644319796441837767922284499546224705478663209861083006852760229926776853038094754552070328034104193400893245372549030298792777754250648573691618715924791001499873681565005082286508541090268441620053266266867900385747199730993494738607288416717003740955088644237161192943830787190330322349014709473604906447330175421589118526590859764383851683929006141834660086562729323700925281924198630248699094612603939433824336106029747514765083", 10);
    
    compute_power_2T_wave(x, T, N1, N2, &saves, out);
}

void multiplication_wave_trial_while (void){
    mpz_t x;
    mpz_t N1;
    mpz_t N2;
    mpz_t out;
    mpz_inits(x, N1, N2, out, NULL);
    vector saves;
    construct_vector(&saves);
    
    // CHANGE PARAMETER HERE
    mpz_set_ui(x, 7);
    mpz_set_str(N1, "48112959837082048697", 10);
    mpz_set_str(N2, "29497513910652490397", 10);
    
    compute_power_2T_wave_while(x, N1, N2, &saves, out);
}


