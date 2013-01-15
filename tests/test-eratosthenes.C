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
#include "../src/smooth-number.h"

using namespace std;

#define NB_LINEAR_RELATIONS	3
#define SIEVING_STEP		10000	//Sieve for SIEVING_STEP numbers at a time

int main(int argc, char **argv)
{
	int modulus_size_bits = 60;
	if(argc >= 2)
	{
		modulus_size_bits = atoi((const char *)argv[1]);
	}

	cout << "[[ Using RSA modulus size " << modulus_size_bits * 2 << "]]" << endl << endl;

	cout << "[Size of int " << sizeof(int) << "]" << endl;
	cout << "[Size of unsgined int " << sizeof(unsigned int) << "]" << endl;
	cout << "[Size of unsgined long int " << sizeof(unsigned long int) << "]" << endl;
	cout << "[Size of unsgined long long int " << sizeof(unsigned long long int) << "]" << endl;
	mpz_class p, q;
	MathUtils::GetRandomPrime(p, modulus_size_bits);
	MathUtils::GetRandomPrime(q, modulus_size_bits);

	cout << "P = " << p << endl;
	cout << "Q = " << q << endl;

	mpz_class N = p * q;
	cout << "N = p*q = " << N << endl << endl;

	mpz_class x = sqrt(N);	//x goes though sqrt(N), sqrt(N)+1, sqrt(N)+2..
	if(x * x == N)			//N is a perfect square, factored!
	{
		cout << "Factored! N = " << x << "²" << endl;
		exit(0);
	}

	SmoothBase smooth_base (N);
	cout << "Smooth Base = " << smooth_base.GetBase() << "\t" << endl;

	TIMER_DECLARE(smooth_base_tm);
	TIMER_START(smooth_base_tm);
		smooth_base.Setup();
	TIMER_STOP(smooth_base_tm);
	cout << "[Number of primes in base " << smooth_base.primes.size() << "]\t";

	TIMER_REPORT(smooth_base_tm);


	vector<SmoothNumber> final_smooth_numbers;

	SmoothNumber *sieving_temp_smooth_numbers;		//We use C arrays here to avoid the copying
													//of pointers (mpz_t) in case of std::vector!
	sieving_temp_smooth_numbers = new SmoothNumber [SIEVING_STEP];

	unsigned long int smooth_base_size = smooth_base.primes.size();
	unsigned long int nb_required_smooth_numbers = smooth_base.primes.size() + NB_LINEAR_RELATIONS;
	unsigned long int nb_discovered_smooth_numbers = 0;

	mpz_class tmp_prime, starting_x_current_round;

	//x goes though sqrt(N), sqrt(N)+1, sqrt(N)+2..
	while(nb_discovered_smooth_numbers < nb_required_smooth_numbers)
	{
		cout << "round " << nb_discovered_smooth_numbers / 20 << endl;

		starting_x_current_round = x;
		//initialize (x, x², exponent vector) for the SIEVING_STEP next x.
		for(int i=0; i<SIEVING_STEP; ++i)
		{
			sieving_temp_smooth_numbers[i].Init(++x, smooth_base_size);
		}

		nb_discovered_smooth_numbers += 20;


		//Reduce sieving_temp_smooth_numbers[] with the primes in the smooth base
		for(unsigned long int i=0; i<smooth_base_size; ++i)
		{
			tmp_prime = smooth_base.primes[i];
			x = starting_x_current_round;
		}
	}

	delete [] sieving_temp_smooth_numbers;

	return 0;
}
