/*
 * eratosthenes.h
 *
 *  Created on: Jan 5, 2013
 *      Author: martani
 *      Copyright martani 2013
 */


#ifndef ERATOSTHENES_H_
#define ERATOSTHENES_H_


#include <vector>
#include <gmpxx.h>
#include <stdint.h>

using namespace std;

class Erastosthenes {
private:

	//For bit at a position x in _primes_bitset: true => x is prime, false => x is not prime
	vector<bool> _primes_bitset;

	//Performs the Eratosthenes sieving for numbers in [0..upperBase]
	//The sieve is saved in the vector _primes_bitset
	void SievePrimesUpTo(unsigned long int upperBase);

	bool _sieving_performed;

	//Indicates whether to free the memory of _primes_bitset immediately after
	//the primes are returned or keep it for later calls
	bool _keep_primes_bitset;

	Erastosthenes(const Erastosthenes& other);
	Erastosthenes& operator= (const Erastosthenes& other);

public:

	Erastosthenes(bool keepMemory = false) : _keep_primes_bitset (keepMemory),
											 _sieving_performed (false)
	{

	}

	~Erastosthenes()
	{
		vector<bool> tmp;
		this->_primes_bitset.clear();
		this->_primes_bitset.swap(tmp);
	}

	//Fills the vector primesArray with the primes in the interval [0..upperBase]
	void GetPrimes(vector<unsigned long int>& primesArray, unsigned long int upperBase);

	//Fills the vector primesArray with the primes in the interval [0..upperBase]
	//and to which N is a quadratic residue
	void GetPrimes_QuadraticResidue(vector<unsigned long int>& primesArray,
			unsigned long int upperBase,
			mpz_class N);
};


#include "eratosthenes.C"

#endif	//ERATOSTHENES_H_
