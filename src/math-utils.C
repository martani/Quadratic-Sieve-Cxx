/*
 * math-utils.C
 *
 *  Created on: Jan 6, 2013
 *      Author: sky
 *      Copyright martani 2013
 */



#ifndef MATH_UTILS_C_
#define MATH_UTILS_C_

#include <sys/time.h>
#include "gmp.h"
#include "mpfr.h"

void MathUtils::GetSmoothnessBase(mpz_class& ret_base, mpz_class& N)
{
	mpfr_t f_N, log_N, log_log_N;
	mpz_t base_mpz;
	mpz_init(base_mpz);

	mpfr_init(f_N); mpfr_init(log_N); mpfr_init(log_log_N);

	mpfr_set_z(f_N, N.get_mpz_t(), MPFR_RNDU);		//f_N = N
	mpfr_log(log_N, f_N, MPFR_RNDU); 				//log_N = log(N)
	mpfr_log(log_log_N, log_N, MPFR_RNDU); 			//log_log_N = log(log(N))

	mpfr_mul(f_N, log_N, log_log_N, MPFR_RNDU); 	//f_N = log(N) * log(log(N))
	mpfr_sqrt(f_N, f_N, MPFR_RNDU); 				//f_N = sqrt(f_N)

	mpfr_div_ui(f_N, f_N, 2, MPFR_RNDU);  			//f_N = f_N/2
	mpfr_exp(f_N, f_N, MPFR_RNDU);					//f_N = e^f_N

	mpfr_get_z(base_mpz, f_N, MPFR_RNDU);
	ret_base = mpz_class(base_mpz);

	mpfr_clears(f_N, log_N, log_log_N, NULL);
}

unsigned long long int MathUtils::GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(unsigned long long int)1000000+tv.tv_usec;
}


void MathUtils::GetRandomPrime(mpz_class& rand_prime, unsigned int num_bits)
{
	gmp_randclass r1 (gmp_randinit_default);
	r1.seed(GetTimeStamp());

	mpz_class rand = r1.get_z_bits(num_bits);

	mpz_t p; mpz_init(p);
	mpz_nextprime(p, rand.get_mpz_t());
	rand_prime = mpz_class(p);

	mpz_clear(p);
}

#endif //MATH_UTILS_C_
