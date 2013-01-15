/*
 * smooth-number.h
 *
 *  Created on: Jan 13, 2013
 *      Author: sky
 */

#ifndef SMOOTH_NUMBER_H_
#define SMOOTH_NUMBER_H_

#include "gmp.h"
#include "gmpxx.h"

class SmoothNumber {
private:
	mpz_t exponent_vector;
	mpz_t x_squared;

	bool _is_initialized;
	unsigned long int smooth_base_size;

	SmoothNumber& operator=(const SmoothNumber& other);
	SmoothNumber (const SmoothNumber& other);

public:
	mpz_class X;

	SmoothNumber ()
	{
		mpz_init(this->x_squared);
		mpz_init(this->exponent_vector);

		this->_is_initialized = true;
		this->smooth_base_size = 0;
	}

	~SmoothNumber ()
	{
		mpz_clear(this->x_squared);
		mpz_clear(this->exponent_vector);
	}


	//Initialize the smooth number with the value of init_x
	//The exponent vector of this smooth number have nb_smooth_primes bits
	//one bit for each prime in the smooth base
	void Init (mpz_class init_x, unsigned long int nb_smooth_primes);


	//Initialize the smooth number without the exponent vector (no smooth base
	//information is saved)
	void InitWithoutExponentVector (mpz_class init_x);

	//Removes all powers of f from x_squared. Same logic as mpz_remove.
	unsigned long int RemovePowerOfFactor(mpz_class f);

	//Returns true if x_squared is a smooth number (=1 after calls to RemovePowerOfFactor)
	bool IsFullyFactoredOnSmoothBase();

	//Sets to 1 the bit corresponding to the i'th prime in the smooth base.
	void SetExponentVectorBit(unsigned long i);
};

void SmoothNumber::Init (mpz_class init_x, unsigned long int nb_smooth_primes)
{
	if(!this->_is_initialized)
	{
		mpz_init(this->x_squared);
		mpz_init2(this->exponent_vector, nb_smooth_primes);

		this->_is_initialized = true;
	}
	else
	{
		//mpz_set(this->x_squared, 0);
		mpz_clear(this->exponent_vector);		//allocate new space for the exponent vector
		mpz_init2(this->exponent_vector, nb_smooth_primes);
	}

	this->X = init_x;
	mpz_mul(this->x_squared, X.get_mpz_t(), X.get_mpz_t());

	this->smooth_base_size = nb_smooth_primes;
}

void SmoothNumber::InitWithoutExponentVector (mpz_class init_x)
{
	this->X = init_x;
	mpz_init(this->x_squared);
	mpz_mul(this->x_squared, X.get_mpz_t(), X.get_mpz_t());

	this->smooth_base_size = 0;
}


unsigned long int SmoothNumber::RemovePowerOfFactor(mpz_class f)
{
	return mpz_remove(this->x_squared, this->x_squared, f.get_mpz_t());
}


bool SmoothNumber::IsFullyFactoredOnSmoothBase()
{
	return mpz_cmp_ui(this->x_squared, 1) == 0;
}


void SmoothNumber::SetExponentVectorBit(unsigned long i)
{
	if(this->smooth_base_size == 0)
			throw std::logic_error("Smooth number must be initialized with the size of the smooth base");

	if(i >= this->smooth_base_size)
		throw std::logic_error("i must be less than the number of primes in the smooth base");

	mpz_setbit(this->exponent_vector, i);
}

#endif /* SMOOTH_NUMBER_H_ */
