/*
 * eratosthenes.C
 *
 *  Created on: Jan 5, 2013
 *      Author: martani
 *      Copyright martani 2013
 */


#ifndef ERATOSTHENES_C_
#define ERATOSTHENES_C_

#include "utils.h"


void Erastosthenes::SievePrimesUpTo(uint64_t upperBase)
{
	//int64_t num_discovered_primes = 0;

	//set all integers from [0..upperBase] as primes in the beginning
	this->_primes_bitset.clear();
	this->_primes_bitset.resize(upperBase + 1, true);

	//mark 0 and 1 as not primes
	this->_primes_bitset[0] = false;
	this->_primes_bitset[1] = false;

	int64_t p = 2, offset;

	//2 is a prime
	//num_discovered_primes = 1;

	while (1) {
			offset = 2 * p;

			//go through all the multiples of the prime p, mark them as composite
			while (offset <= upperBase) {
				this->_primes_bitset[offset] = false;
				offset += p;
			}

			//move to next prime
			offset = p + 1;

			//search for the next next prime (ie. non marked integer in _primes_bitset)
			while (offset <= upperBase && !this->_primes_bitset[offset]) {
				offset++;
			}

			if (offset > upperBase)
				break;

			//next prime is in offset now
			p = offset;
			//num_discovered_primes++;
	}

	this->_sieving_performed = true;
}

void Erastosthenes::GetPrimes(vector<uint64_t>& primesArray, uint64_t upperBase)
{
	vector<uint64_t> primes;

	//Sieving
	if(!this->_sieving_performed)
		Erastosthenes::SievePrimesUpTo(upperBase);

	//For bit at a position x in _primes_bitset: true => x is prime, false => x is not prime
	for(uint64_t i=0; i<this->_primes_bitset.size(); ++i)
	{
		if(this->_primes_bitset[i])
			primes.push_back(i);
	}

	primesArray.swap(primes);

	if(!this->_keep_primes_bitset)
	{
		vector<bool> tmp;
		this->_primes_bitset.clear();
		this->_primes_bitset.swap(tmp);

		this->_sieving_performed = false;
	}
}

void Erastosthenes::GetPrimes_QuadraticResidue(vector<uint64_t>& primesArray,
		uint64_t upperBase, mpz_class N)
{
	vector<uint64_t> primes;

	//Get an mpz_t out of C++ mpz_class
	mpz_t Nmpz;
	mpz_init_set(Nmpz, N.get_mpz_t());

	//p is an mpz_t to hold all the primes
	mpz_t p; mpz_init(p);

	//Sieving
	if(!this->_sieving_performed)
		Erastosthenes::SievePrimesUpTo(upperBase);

	//Push the prime 2, mpz_legendre() required an odd positive prime
	//primes.push_back(2);

	//For bit at a position x in _primes_bitset: true => x is prime, false => x is not prime
	for(uint64_t i=0; i<this->_primes_bitset.size(); ++i)
	{
		mpz_set_ui(p, i);
		if(this->_primes_bitset[i] && mpz_legendre(Nmpz, p) == 1)
			primes.push_back(i);
	}

	primesArray.swap(primes);

	if(!this->_keep_primes_bitset)
	{
		vector<bool> tmp;
		this->_primes_bitset.clear();
		this->_primes_bitset.swap(tmp);

		this->_sieving_performed = false;
	}
}

#endif //ERATOSTHENES_C_
