/*
 * QS-parallel.h
 *
 *  Created on: Mar 5, 2013
 *      Author: sky
 */

#ifndef QS_PARALLEL_H_
#define QS_PARALLEL_H_

#include <omp.h>
#include "QS.h"

class QSParallel : public QS {

private:
	unsigned nb_threads;

	virtual void Sieve ();


public:

	QSParallel (mpz_class N) : QS (N)
	{
		this->nb_threads = omp_get_max_threads ();
	}

	QSParallel (mpz_class N, unsigned nb_additional_linear_relations)
	: QS (N, nb_additional_linear_relations)
	{
		this->nb_threads = omp_get_max_threads ();
	}

	/**
	 * Sets the number of threads to use in the sieving step.
	 * If this method is not called, the max available number of threads
	 * on this processor are used.
	 */
	void setNumThreads (unsigned nbThreads);
	virtual ~QSParallel () { ; }
};


void QSParallel::setNumThreads (unsigned nbThreads)
{
	if(nbThreads < 1)
		throw logic_error ("nbThreads must be at least 1");

	this->nb_threads = nbThreads;
}


void QSParallel::Sieve()
{
	//Set num threads.
	omp_set_num_threads(this->nb_threads);
	dcout << "[[Parallel Sieve() using " << this->nb_threads << " threads.]]\n";

	mpz_class x = sqrt(this->_N);	//x goes though sqrt(N), sqrt(N)+1, sqrt(N)+2..
	if(x * x == this->_N)			//N is a perfect square, factored!
	{
		dcout << "Factored! N = " << x << "²" << "\n";
		SetFactors(x, x);
		return;
	}
	++x;		//Start from sqrt(N)+1

	SmoothBase smooth_base (this->_N);
	dcout << "Smooth Base = " << smooth_base.GetBase() << "\t\n";

	this->SetupSmoothBase(smooth_base);

	this->smooth_base_size = smooth_base.primes.size ();
	unsigned long smooth_base_size = smooth_base.primes.size ();
	unsigned long nb_required_smooth_numbers = smooth_base.primes.size () + this->additional_linear_relations;
	unsigned long nb_discovered_smooth_numbers = 0;

	//The exponents vectors' matrix
	this->_M.Init(nb_required_smooth_numbers, smooth_base_size);

	std::vector<SmoothNumber> sieving_temp_smooth_numbers (SIEVING_STEP);

	mpz_t tmp_prime_p; mpz_init(tmp_prime_p);  //using mpz_t for performance
	unsigned long exponent_prime_p;
	mpz_class starting_x_current_round;  //goes through sqrt(N), sqrt(N)+SIEVING_STEP ...
	mpz_class x_idx;

	TIMER_DECLARE(sieving_timer);
	TIMER_START(sieving_timer);
	dcout << "\n";
	dcout << "\n";

	starting_x_current_round = x;

	//x goes though sqrt(N), sqrt(N)+1, sqrt(N)+2..
	while(nb_discovered_smooth_numbers < nb_required_smooth_numbers)
	{
		x = starting_x_current_round;

		dcout << "\r\tSieving at " << starting_x_current_round << "\t";
		dcout << "Smooth numbers found\t " << nb_discovered_smooth_numbers
			  << "/" << smooth_base_size << "          " << std::ends;
		dcout.flush();

		//initialize (x, x², exponent vector) for the SIEVING_STEP next x.
		for(int i=0; i<SIEVING_STEP; ++i)
		{
			sieving_temp_smooth_numbers[i].Init(x++, smooth_base_size, this->_N);
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

#pragma omp parallel for private (exponent_prime_p)
			for(unsigned long j=x_idx.get_ui(); j<SIEVING_STEP;
					j+=smooth_base.primes[i])
			{
				//Eliminate all powers of tmp_prime in sieving_temp_smooth_numbers[j]
				exponent_prime_p = sieving_temp_smooth_numbers[j].RemovePowerOfFactor(tmp_prime_p);

				//sieving_temp_smooth_numbers[j] is guaranteed to be divisible by tmp_prime_p
				//TODO: Comment in release code
				assert(exponent_prime_p != 0);

				//If the power of tmp_prime_p is odd, we set the corresponding bit in the exponent vector
				if(exponent_prime_p & 1)
				{
					sieving_temp_smooth_numbers[j].SetExponentVectorBit(i);
				}

				/*if(exponent_prime_p != 0)
				{
					sieving_temp_smooth_numbers[j].SetNonZeroExponentsVectorBit(i);
				}*/

				//Can never happen
				/*if(exponent_prime_p == 0 &&
						sieving_temp_smooth_numbers[j].IsFullyFactoredOnSmoothBase())
					cout << " REPETITION " << endl;*/

				//Check if sieving_temp_smooth_numbers[j] was fully factored
				if(sieving_temp_smooth_numbers[j].IsFullyFactoredOnSmoothBase())
				{
					#pragma omp atomic
					++nb_discovered_smooth_numbers;

					//Add the exponent vector to the matrix
					this->_M.PushExponentVector(sieving_temp_smooth_numbers[j]);

					//MATRIX_STATS.PushRow(sieving_temp_smooth_numbers[j].GetNonZeroExponentsVector ());

					//Add the smooth number to the list
					SmoothNumber tmp;
					final_smooth_numbers.push_back(tmp);
					final_smooth_numbers.back ()
										.InitWithoutExponentVector(sieving_temp_smooth_numbers[j].X, this->_N);
				}
			}

/////////// Second root of (N mod (smooth_base.primes[i]))
			if(smooth_base.roots_1[i] == smooth_base.roots_2[i])
				continue;

			x = MathUtils::GetNextMultipleGreaterThanX(starting_x_current_round,
								mpz_class(tmp_prime_p), smooth_base.roots_2[i]);

			x_idx = x-starting_x_current_round;

#pragma omp parallel for private (exponent_prime_p)
			for(unsigned long int j= x_idx.get_ui(); j<SIEVING_STEP; j+=smooth_base.primes[i])
			{
				//Eliminate all powers of tmp_prime in sieving_temp_smooth_numbers[j]
				exponent_prime_p = sieving_temp_smooth_numbers[j].RemovePowerOfFactor(tmp_prime_p);

				//sieving_temp_smooth_numbers[j] is guaranteed to be divisible by tmp_prime_p
				//TODO: Comment in release code
				assert(exponent_prime_p != 0);

				//If the power of tmp_prime_p is odd, we set the corresponding bit in the exponent vector
				if(exponent_prime_p & 1)
				{
					sieving_temp_smooth_numbers[j].SetExponentVectorBit(i);
				}

				/*if(exponent_prime_p != 0)
				{
					sieving_temp_smooth_numbers[j].SetNonZeroExponentsVectorBit(i);
				}*/

				//Can never happen
				/*if(exponent_prime_p == 0 &&
						sieving_temp_smooth_numbers[j].IsFullyFactoredOnSmoothBase())
					cout << " ERROR " << endl;*/

				//Check if sieving_temp_smooth_numbers[j] was fully factored
				if(sieving_temp_smooth_numbers[j].IsFullyFactoredOnSmoothBase())
				{
					#pragma omp atomic
					++nb_discovered_smooth_numbers;

					//Add the exponent vector to the matrix
					this->_M.PushExponentVector(sieving_temp_smooth_numbers[j]);

					//MATRIX_STATS.PushRow(sieving_temp_smooth_numbers[j].GetNonZeroExponentsVector ());

					//Add the smooth number to the list
					//Add the smooth number to the list
					final_smooth_numbers.push_back(SmoothNumber ());
					final_smooth_numbers.back ()
										.InitWithoutExponentVector(sieving_temp_smooth_numbers[j].X, this->_N);
				}
			}

		}  //All primes in the smooth base are handled at this point

		//Move to the next interval
		starting_x_current_round += SIEVING_STEP;
	}

	mpz_clear(tmp_prime_p);
	dcout << "\r\tSieving at " << starting_x_current_round << "\t";
	dcout << "Smooth numbers found\t " << nb_discovered_smooth_numbers
		  << "/" << smooth_base_size << "          " << ends;
	cout.flush();
	dcout << "\n";
	dcout << "\n";

	TIMER_STOP(sieving_timer);
	TIMER_REPORT(sieving_timer);

	//Utils::dumpMatrixAsPbmImage(MATRIX_STATS, "stats.pbm");
}



#endif /* QS_PARALLEL_H_ */
