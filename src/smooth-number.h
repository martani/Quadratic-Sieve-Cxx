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
public:
	mpz_class X;

	SmoothNumber ()
	{
		mpz_init(this->x_squared);
		mpz_init(this->exponent_vector);
//		mpz_init(this->non_zeor_exponents_vector);

		this->_is_initialized = true;
		this->smooth_base_size = 0;
	}

	SmoothNumber (const SmoothNumber& other)
	{
		this->X = other.X;
		mpz_init_set(this->x_squared, other.x_squared);
		mpz_init_set(this->exponent_vector, other.exponent_vector);
//		mpz_init(this->non_zeor_exponents_vector);

		this->_is_initialized = other._is_initialized;
		this->smooth_base_size = other.smooth_base_size;
	}

	SmoothNumber& operator=(const SmoothNumber& other)
	{
		this->X = other.X;
		mpz_init_set(this->x_squared, other.x_squared);
		mpz_init_set(this->exponent_vector, other.exponent_vector);
//		mpz_init(this->non_zeor_exponents_vector);

		this->_is_initialized = other._is_initialized;
		this->smooth_base_size = other.smooth_base_size;

		return *this;
	}

	~SmoothNumber ()
	{
		mpz_clear(this->x_squared);
		mpz_clear(this->exponent_vector);
//		mpz_clear(this->non_zeor_exponents_vector);
	}


	//Initialize the smooth number with the value of init_x
	//The exponent vector of this smooth number have nb_smooth_primes bits
	//one bit for each prime in the smooth base
	inline void Init (mpz_class init_x, unsigned long int nb_smooth_primes, mpz_class modulus);


	//Initialize the smooth number without the exponent vector (no smooth base
	//information is saved)
	void InitWithoutExponentVector (mpz_class init_x, mpz_class modulus);

	//Removes all powers of f from x_squared. Same logic as mpz_remove.
	inline unsigned long int RemovePowerOfFactor(mpz_t f);

	//Returns true if x_squared is a smooth number (=1 after calls to RemovePowerOfFactor)
	inline bool IsFullyFactoredOnSmoothBase();

	//Sets to 1 the bit corresponding to the i'th prime in the smooth base.
	inline void SetExponentVectorBit(unsigned long i);

	inline void SetNonZeroExponentsVectorBit(unsigned long i);

	const mpz_t& GetExponentVector ();
//	const mpz_t& GetNonZeroExponentsVector ();

	mpz_class GetXSquared () const;


private:
	mpz_class modulus_N;
	mpz_t x_squared;
	mpz_t exponent_vector;
//	mpz_t non_zeor_exponents_vector;	//bit is set to one if corresponding prime is <> 0

	bool _is_initialized;
	unsigned long int smooth_base_size;
};

inline void SmoothNumber::Init (mpz_class init_x, unsigned long int nb_smooth_primes,
		mpz_class modulus)
{
	if(!this->_is_initialized)
	{
		mpz_init(this->x_squared);
		mpz_init2(this->exponent_vector, nb_smooth_primes);
		mpz_set_ui(this->exponent_vector, 0);

//		mpz_init2(this->non_zeor_exponents_vector, nb_smooth_primes);
//		mpz_set_ui(this->non_zeor_exponents_vector, 0);

		this->_is_initialized = true;
	}
	else
	{
		//mpz_set(this->x_squared, 0);
		mpz_clear(this->exponent_vector);		//allocate new space for the exponent vector
		mpz_init2(this->exponent_vector, nb_smooth_primes);
		mpz_set_ui(this->exponent_vector, 0);

//		mpz_clear(this->non_zeor_exponents_vector);		//allocate new space for the exponent vector
//		mpz_init2(this->non_zeor_exponents_vector, nb_smooth_primes);
//		mpz_set_ui(this->non_zeor_exponents_vector, 0);
	}

	this->X = init_x;
	//Set x_squared to X*X - N
	mpz_pow_ui(this->x_squared, this->X.get_mpz_t (), 2);
	mpz_sub(this->x_squared, this->x_squared, modulus.get_mpz_t ());

	this->smooth_base_size = nb_smooth_primes;
	this->modulus_N = modulus;
}

void SmoothNumber::InitWithoutExponentVector (mpz_class init_x, mpz_class modulus)
{
	this->X = init_x;
	if(!this->_is_initialized)
	{
		mpz_init(this->x_squared);
		this->_is_initialized = true;
	}

	//Set x_squared to X*X - N
	mpz_pow_ui(this->x_squared, X.get_mpz_t (), 2);
	mpz_sub(this->x_squared, this->x_squared, modulus.get_mpz_t ());

	this->smooth_base_size = 0;
	this->modulus_N = modulus;
}


inline unsigned long int SmoothNumber::RemovePowerOfFactor(mpz_t f)
{
	return mpz_remove(this->x_squared, this->x_squared, f);
}


inline bool SmoothNumber::IsFullyFactoredOnSmoothBase()
{
	return mpz_cmp_ui(this->x_squared, 1) == 0;
}


inline void SmoothNumber::SetExponentVectorBit(unsigned long i)
{
//	if(this->smooth_base_size == 0)
//			throw std::logic_error("Smooth number must be initialized with the size of the smooth base");
//
//	if(i >= this->smooth_base_size)
//		throw std::logic_error("i must be less than the number of primes in the smooth base");

	mpz_setbit(this->exponent_vector, i);
}

//inline void SmoothNumber::SetNonZeroExponentsVectorBit(unsigned long i)
//{
//	mpz_setbit(this->non_zeor_exponents_vector, i);
//}

const mpz_t& SmoothNumber::GetExponentVector ()
{
	return this->exponent_vector;
}

mpz_class SmoothNumber::GetXSquared () const
{
	return mpz_class (this->x_squared);
}
//const mpz_t& SmoothNumber::GetNonZeroExponentsVector () const
//{
//	return this->non_zeor_exponents_vector;
//}

#endif /* SMOOTH_NUMBER_H_ */
