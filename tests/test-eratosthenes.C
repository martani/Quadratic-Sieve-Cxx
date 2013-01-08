/*
 * test-eratosthenes.C
 *
 *  Created on: Jan 5, 2013
 *      Author: martani
 *      Copyright martani 2013
 */


#include <iostream>
#include <cstdlib>
#include "../src/eratosthenes.h"
#include "../src/utils.h"
#include "../src/math-utils.h"
#include "../src/smooth-base.h"

using namespace std;

int main(int argc, char **argv)
{
	int modulus_size_bits = 60;
	if(argc >= 2)
	{
		modulus_size_bits = atoi((const char *)argv[1]);
	}

	cout << "[[ Using RSA modulus size " << modulus_size_bits * 2 << "]]" << endl << endl;

	cout << "[Size of unsgined long int " << sizeof(unsigned long int) << "]" << endl;
	cout << "[Size of unsgined long long int " << sizeof(unsigned long long int) << "]" << endl;
	mpz_class p, q;
	MathUtils::GetRandomPrime(p, modulus_size_bits);
	MathUtils::GetRandomPrime(q, modulus_size_bits);

	cout << "P = " << p << endl;
	cout << "Q = " << q << endl;

	mpz_class N = p*q;
	cout << "N = p*q = " << N << endl << endl;

	TIMER_DECLARE(smooth_base_tm);
	TIMER_START(smooth_base_tm);

		SmoothBase smooth_base (N);
		smooth_base.Setup();

	TIMER_STOP(smooth_base_tm);

	cout << "Smooth Base = " << smooth_base.GetBase() << "\t";
	cout << "[Number of primes in base " << smooth_base.primes.size() << "]\t";

	TIMER_REPORT(smooth_base_tm);

	return 0;
}
