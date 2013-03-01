/*
 * QS.C
 *
 *  Created on: Feb 21, 2013
 *      Author: sky
 */


#ifndef QS_C_
#define QS_C_

#include <vector>
#include <iostream>
#include <cassert>

#include "QS.h"



QS::QS (mpz_class N)
{
	_N = N;
	p = 1;
	q = N;

	this->additional_linear_relations = 5;	//5 is _probably_ enough
	this->smooth_base_size = 0;
	this->computation_completed = false;
}

QS::QS (mpz_class N, unsigned nb_additional_linear_relations)
{
	_N = N;
	p = 1;
	q = N;

	//20 more relations is estimated as sufficient
	if(nb_additional_linear_relations < MAX_NB_ADDITIONAL_LINEAR_RELS)
		this->additional_linear_relations =
				nb_additional_linear_relations == 0 ?
						1 :
						nb_additional_linear_relations;
	else
		this->additional_linear_relations = MAX_NB_ADDITIONAL_LINEAR_RELS;

	this->smooth_base_size = 0;
	this->computation_completed = false;
}


void QS::Factor()
{
	//Step 1: Sieve for smooth numbers
	try
	{
		this->Sieve ();
	}
	catch(std::exception &e)
	{
		std::cout << "[[ERROR]] Cannot setup smooth base, the integer is too large" << endl
				  << "Exception message:\t"
				  << e.what () << endl;

		return;
	}

	//Step 2: Perform linear algebra
	this->PerformeGaussianElimination ();

	//Step 3: Factor using results from linear algebra
	this->FactorUsingLinearRelations ();
}

void QS::Sieve()
{
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

	TIMER_DECLARE(smooth_base_tm);
	TIMER_START(smooth_base_tm);
		smooth_base.Setup();
	TIMER_STOP(smooth_base_tm);
	dcout << "[Number of primes in base " << smooth_base.primes.size() << "]\t";
	TIMER_REPORT(smooth_base_tm);

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


void QS::PerformeGaussianElimination()
{
	dcout << "\n";
	dcout << "<< Performing Linear Algebra >>" << "\n";
	TIMER_DECLARE(gauss);
	TIMER_START(gauss);

	this->gauss.Echelonize(this->_M);
	this->linear_relations = gauss.GetLinearRelations();

	TIMER_STOP(gauss);
	TIMER_REPORT(gauss);

	dcout << "\tLinear dependent relations found : " << this->linear_relations.size ()
		  << "\n";
	//Utils::dumpMatrixAsPbmImage(M, "M.pbm");
}


void QS::FactorUsingLinearRelations()
{
	dcout << "\n";
	dcout << "<< Factoring >>\n";
	mpz_class x_side, y_side, factor;
	mpz_t fctr; mpz_init(fctr);

	for(int i = 0; i < this->gauss.GetNbLinearRelations (); ++i)
	{
		dcout << "Trying relation " << i << "\n";
		x_side = 1;
		y_side = 1;

		//For all bits representing the primes in the smooth base
		for(int j = 0; j < this->smooth_base_size; ++j)
		{
			if (mpz_tstbit(this->linear_relations[i].get_mpz_t (), j))
			{
				x_side *= final_smooth_numbers[j].X;  //j'th smooth number is used in the product
				x_side %= this->_N;  //Reduce modulo N to keep the numbers fairly small

				y_side *= final_smooth_numbers[j].GetXSquared ();  //We cannot reduce that!!
			}
		}

		//cout << "Y SIDE IS " << y_side <<  endl;

		y_side = sqrt(y_side);  //
		y_side = y_side % this->_N;  //Now we can reduce by N

		x_side = x_side - y_side;


		mpz_gcd(fctr, x_side.get_mpz_t (), this->_N.get_mpz_t ());

		factor = mpz_class (fctr);
		if(factor != this->_N && factor != 1)  //If the factor is not trivial, N or 1, then we got it!
			break;
	}
	mpz_clear(fctr);

	if(factor == 1 || factor == this->_N)
	{
		dcout << ">>>> Failed to factor " << this->_N << " <<<<\t"
				<< "Try using more linear relations" << "\n";
	}
	else
	{
		dcout << "\n";
		dcout << ">>>>>>> Factored " << this->_N << "\n";
		dcout << "\t Factor 1: " << factor << "\n";
		dcout << "\t Factor 2: " << this->_N/factor << "\n";

		SetFactors(factor, this->_N/factor);
	}
}

void QS::SetFactors(mpz_class f1, mpz_class f2)
{
	this->p = f1;
	this->q = f2;

	this->computation_completed = true;
}

mpz_class QS::GetFactor1()
{
	return this->p;
}

mpz_class QS::GetFactor2()
{
	return this->q;
}

#endif /* QS_C_ */

