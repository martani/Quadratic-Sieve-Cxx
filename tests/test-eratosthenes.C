/*
 * test-eratosthenes.C
 *
 *  Created on: Jan 5, 2013
 *      Author: martani
 *      Copyright martani 2013
 */


#include <iostream>
#include <cstdlib>
#include <cassert>
#include "../src/eratosthenes.h"
#include "../src/utils.h"
#include "../src/math-utils.h"
#include "../src/smooth-base.h"
#include "../src/smooth-number.h"
#include "../src/linear-algebra/matrix.h"
#include "../src/linear-algebra/gauss-elimination.h"

#include "../src/debug-cout.h"

using namespace std;

#define NB_LINEAR_RELATIONS	10
#define SIEVING_STEP		50000	//Sieve for SIEVING_STEP numbers at a time

int main(int argc, char **argv)
{
	int modulus_size_bits = 50;
	if(argc >= 2)
	{
		//modulus_size_bits = atoi((const char *)argv[1]);
	}

	cout << "[[GMP " << gmp_version << "]"
		 << " [MPFR " << mpfr_get_version() << "]]" << endl;
	cout << "[[ Using RSA modulus size " << modulus_size_bits * 2 << "]]" << endl << endl;

	cout << "[Size of int " << sizeof(int) << "]" << endl;
	cout << "[Size of unsgined int " << sizeof(unsigned int) << "]" << endl;
	cout << "[Size of unsgined long int " << sizeof(unsigned long int) << "]" << endl;
	cout << "[Size of unsgined long long int " << sizeof(unsigned long long int) << "]" << endl << endl;
	mpz_class p, q;
	MathUtils::GetRandomPrime(p, modulus_size_bits);
	MathUtils::GetRandomPrime(q, modulus_size_bits);

	mpz_class N;
	if(argc >= 2)
	{
		//modulus_size_bits = atoi((const char *)argv[1]);
		N = argv[1];
	}
	else
	{
		cout << "P = " << p << endl;
		cout << "Q = " << q << endl;

		N = p * q;
	}
	cout << "N = P*Q = " << N << endl << endl;


	mpz_class x = sqrt(N);	//x goes though sqrt(N), sqrt(N)+1, sqrt(N)+2..
	if(x * x == N)			//N is a perfect square, factored!
	{
		cout << "Factored! N = " << x << "²" << endl;
		exit(0);
	}
	++x;		//Start from sqrt(N)+1

	SmoothBase smooth_base (N);
	cout << "Smooth Base = " << smooth_base.GetBase() << "\t" << endl;

	TIMER_DECLARE(smooth_base_tm);
	TIMER_START(smooth_base_tm);
		smooth_base.Setup();
	TIMER_STOP(smooth_base_tm);
	cout << "[Number of primes in base " << smooth_base.primes.size() << "]\t";
	TIMER_REPORT(smooth_base_tm);

//	Utils::AppendVectorToFile(smooth_base.primes, "primes_roots.txt");
//	Utils::AppendVectorToFile(smooth_base.roots_1, "primes_roots.txt");
//	Utils::AppendVectorToFile(smooth_base.roots_2, "primes_roots.txt");
//	assert (smooth_base.primes.front() == 2);

	unsigned long int smooth_base_size = smooth_base.primes.size();
	unsigned long int nb_required_smooth_numbers = smooth_base.primes.size() + NB_LINEAR_RELATIONS;
	unsigned long int nb_discovered_smooth_numbers = 0;


	Matrix M (nb_required_smooth_numbers, smooth_base_size);
//	Matrix MATRIX_STATS (nb_required_smooth_numbers, smooth_base_size);

	vector<SmoothNumber> final_smooth_numbers;
	vector<SmoothNumber> sieving_temp_smooth_numbers (SIEVING_STEP);		//We use C arrays here to avoid the copying
													//of pointers (mpz_t) in case of std::vector!
	//sieving_temp_smooth_numbers = new SmoothNumber [SIEVING_STEP];

	mpz_t tmp_prime_p; mpz_init(tmp_prime_p);			//using mpz_t for performance
	unsigned long int exponent_prime_p;
	mpz_class starting_x_current_round;		//goes throught sqrt(N), sqrt(N)+SIEVING_STEP ...
	mpz_class x_idx;

	TIMER_DECLARE(sieving_timer);
	TIMER_START(sieving_timer);
	cout << endl << endl;
	starting_x_current_round = x;


	//x goes though sqrt(N), sqrt(N)+1, sqrt(N)+2..
	while(nb_discovered_smooth_numbers < nb_required_smooth_numbers)
	{
		x = starting_x_current_round;

		cout << "\r\tSieving at " << starting_x_current_round << "\t";
		cout << "Smooth numbers found\t " << nb_discovered_smooth_numbers
			 << "/" << smooth_base_size << "          " << ends;
		cout.flush();

		//initialize (x, x², exponent vector) for the SIEVING_STEP next x.
		for(int i=0; i<SIEVING_STEP; ++i)
		{
			sieving_temp_smooth_numbers[i].Init(x++, smooth_base_size, N);
		}

		//Reduce sieving_temp_smooth_numbers[] with the primes in the smooth base
		for(unsigned long int i=0; i<smooth_base_size; ++i)
		{
			//tmp_prime_p = smooth_base.primes[i];
			mpz_set_ui(tmp_prime_p, smooth_base.primes[i]);

/////////// First root of (N mod (smooth_base.primes[i]))
			x = MathUtils::GetNextMultipleGreaterThanX(starting_x_current_round,
					mpz_class(tmp_prime_p), smooth_base.roots_1[i]);

			//assert(x >= starting_x_current_round);

			x_idx = x-starting_x_current_round;

			for(unsigned long int j= x_idx.get_ui(); j<SIEVING_STEP; j+=smooth_base.primes[i])
			{
				//Eliminate all powers of tmp_prime in sieving_temp_smooth_numbers[j]
				exponent_prime_p = sieving_temp_smooth_numbers[j].RemovePowerOfFactor(tmp_prime_p);
				assert(exponent_prime_p != 0);

				//If the power of tmp_prime_p is odd, we set the corresponding bit in the exponent vector
				if(exponent_prime_p & 1)
				{
					sieving_temp_smooth_numbers[j].SetExponentVectorBit(i);
				}

//				if(exponent_prime_p != 0)
//				{
//					sieving_temp_smooth_numbers[j].SetNonZeroExponentsVectorBit(i);
//				}

				//Can never happen
//				if(exponent_prime_p == 0 && sieving_temp_smooth_numbers[j].IsFullyFactoredOnSmoothBase())
//					cout << " REPETITION " << endl;

				//Check if sieving_temp_smooth_numbers[j] was fully factored
				if(sieving_temp_smooth_numbers[j].IsFullyFactoredOnSmoothBase())
				{
					++nb_discovered_smooth_numbers;

					//Add the exponent vector to the matrix
					M.PushExponentVector(sieving_temp_smooth_numbers[j]);
//					MATRIX_STATS.PushRow(sieving_temp_smooth_numbers[j].GetNonZeroExponentsVector ());
					//Add the smooth number to the list
					final_smooth_numbers.push_back(SmoothNumber ());
					final_smooth_numbers.back ().InitWithoutExponentVector(sieving_temp_smooth_numbers[j].X, N);
				}
			}


/////////// Second root of (N mod (smooth_base.primes[i]))
			if(smooth_base.roots_1[i] == smooth_base.roots_2[i])
				continue;

			x = MathUtils::GetNextMultipleGreaterThanX(starting_x_current_round,
					mpz_class(tmp_prime_p), smooth_base.roots_2[i]);

			x_idx = x-starting_x_current_round;
			for(unsigned long int j= x_idx.get_ui(); j<SIEVING_STEP; j+=smooth_base.primes[i])
			{
				//Eliminate all powers of tmp_prime in sieving_temp_smooth_numbers[j]
				exponent_prime_p = sieving_temp_smooth_numbers[j].RemovePowerOfFactor(tmp_prime_p);
				assert(exponent_prime_p != 0);
				//If the power of tmp_prime_p is odd, we set the corresponding bit in the exponent vector
				if(exponent_prime_p & 1)
				{
					sieving_temp_smooth_numbers[j].SetExponentVectorBit(i);
				}

//				if(exponent_prime_p != 0)
//				{
//					sieving_temp_smooth_numbers[j].SetNonZeroExponentsVectorBit(i);
//				}

//				if(exponent_prime_p == 0 && sieving_temp_smooth_numbers[j].IsFullyFactoredOnSmoothBase())
//					cout << " REPETITION " << endl;

				//Check if sieving_temp_smooth_numbers[j] was fully factored
				if(sieving_temp_smooth_numbers[j].IsFullyFactoredOnSmoothBase())
				{
					++nb_discovered_smooth_numbers;

					//Add the exponent vector to the matrix
					M.PushExponentVector(sieving_temp_smooth_numbers[j]);
//					MATRIX_STATS.PushRow(sieving_temp_smooth_numbers[j].GetNonZeroExponentsVector ());
					//Add the smooth number to the list
					final_smooth_numbers.push_back(SmoothNumber ());
					final_smooth_numbers.back ().InitWithoutExponentVector(sieving_temp_smooth_numbers[j].X, N);
				}
			}
		}	//handled of all primes

		//Move to the next interval
		starting_x_current_round += SIEVING_STEP;
	}
	cout << "\r\tSieving at " << starting_x_current_round << "\t";
	cout << "Smooth numbers found\t " << nb_discovered_smooth_numbers
		 << "/" << smooth_base_size << "          " << ends;
	cout.flush();

	cout << endl << endl;

	TIMER_STOP(sieving_timer);
	TIMER_REPORT(sieving_timer);

//	Utils::dumpMatrixAsPbmImage(MATRIX_STATS, "stats.pbm");

//////// Linear Algebra
	cout << endl << "<< Performing Linear Algebra >>" << endl;
	TIMER_DECLARE(gauss);
	TIMER_START(gauss);

	GaussElimination gauss;
	gauss.Echelonize(M);
	mpz_t *linear_relations = gauss.GetLinearRelations();

	TIMER_STOP(gauss);
	TIMER_REPORT(gauss);

	Utils::dumpMatrixAsPbmImage(M, "M.pbm");

/////// Factoring
	cout << endl << "<< Factoring >>" << endl;
	mpz_class x_side, y_side, factor;
	for(int i=0; i<gauss.GetNbLinearRelations (); ++i)
	{
		cout << "Trying relation " << i << endl;
		x_side = 1;
		y_side = 1;

		//For all bits representing the primes in the smooth base
		for(int j=0; j<smooth_base_size; ++j)
		{
			if (mpz_tstbit(linear_relations[i], j))
			{
				x_side *= final_smooth_numbers[j].X;  //j'th smooth number is used in the product
				x_side %= N;  //Reduce modulo N to keep the numbers fairly small

				y_side *= final_smooth_numbers[j].GetXSquared ();  //We cannot reduce that!!
			}
		}

		//cout << "Y SIDE IS " << y_side <<  endl;

		y_side = sqrt(y_side);  //
		y_side = y_side % N;  //Now we can reduce by N

		x_side = x_side - y_side;

		mpz_t fctr; mpz_init(fctr);
		mpz_gcd(fctr, x_side.get_mpz_t (), N.get_mpz_t ());

		factor = mpz_class (fctr);
		if(factor != N && factor != 1)  //If the factor is not trivial, N or 1, then we got it!
			break;

		mpz_clear(fctr);
	}

	if(factor == 1 || factor == N)
	{
		cout << ">>>> Failed to factor " << N << " <<<<\t"
				<< "Try using more linear relations" << endl;
	}
	else
	{
		cout << endl << ">>>>>>> Factored " << N << endl;
		cout << "\t Factor 1: " << factor << endl
			 << "\t Factor 2: " << N/factor << endl;
	}

	return 0;
}
