/*
 * math-utils.h
 *
 *  Created on: Jan 6, 2013
 *      Author: sky
 *      Copyright martani 2013
 */


#ifndef MATH_UTILS_H_
#define MATH_UTILS_H_

#include <gmpxx.h>


class MathUtils
{
public:

	//Computes the smoothness base to which all the primes used in sieving
	//are in [0..rest_base]
	static void GetSmoothnessBase(mpz_class& ret_base, mpz_class& n);

	static void GetRandomPrime(mpz_class& rand_prime, unsigned int num_bits);

	static unsigned long int GetTimeStamp();
};


#include "math-utils.C"

#endif //MATH_UTILS_H_
