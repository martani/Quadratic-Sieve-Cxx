/*
 * test-qs.C
 *
 *  Created on: Feb 28, 2013
 *      Author: sky
 */


#include <iostream>
#include <algorithm>
#include <gmp.h>
#include <gmpxx.h>
#include "../src/QS.h"
#include "../src/QS-parallel.h"

using namespace std;

#define NB_LINEAR_RELATIONS	10

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return NULL;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

void printInfo()
{
	cout << "[[GMP " << gmp_version << "]"
		 << " [MPFR " << mpfr_get_version() << "]]" << endl;

	cout << "[Size of int " << sizeof(int) << "]" << endl;
	cout << "[Size of unsgined int " << sizeof(unsigned int) << "]" << endl;
	cout << "[Size of unsgined long int " << sizeof(unsigned long int) << "]" << endl;
	cout << "[Size of unsgined long long int " << sizeof(unsigned long long int) << "]" << endl << endl;
}

void printUsage(string program_name)
{
	cout << "Usage: " << endl
		 << program_name << " -b B"<< "\t" << "-Generates a random RSA number of size 2*B bits"
		 << endl << "\t\t (each prime of size B bits) and factors it."
		 << endl << "\t\t NB: B must be at least 20 bits, it has very poor performance with small integers."
		 << endl << endl
		 << program_name << " N" << "\t" << "-Factors N (N must be an RSA number)."
		 << endl << endl;
}

int main(int argc, char **argv)
{
	//the RSA integer to factor
	mpz_class N = 0;

	//the factors
	mpz_class p, q;

	//Parallel sieving?
	bool parallel = false;
	int nb_parallel_threads = 0;

	if(argc < 2 || cmdOptionExists(argv, argv+argc, "-h"))
	{
		printUsage(argv[0]);
		return 0;
	}

	//generate a random RSA integer?
	int modulus_size_bits;
	if(cmdOptionExists(argv, argv+argc, "-b"))
	{
		char *prime_size = getCmdOption(argv, argv+argc, "-b");
		if(prime_size == NULL)
		{
			printUsage(argv[0]);
			return -1;
		}

		modulus_size_bits = atoi(prime_size);
		if(modulus_size_bits < 20)
		{
			printUsage(argv[0]);
			return -1;
		}

		cout << "Generating random primes of size " << modulus_size_bits
			 << "..." << endl;
		MathUtils::GetRandomPrime(p, modulus_size_bits);
		MathUtils::GetRandomPrime(q, modulus_size_bits);

		N = p*q;
		cout << "P = " << p << endl;
		cout << "Q = " << q << endl;
	}

	//is the RSA integer provided?
	if(!cmdOptionExists(argv, argv+argc, "-b") && argc >= 2)
	{
		N = argv[1];
	}

	if(N == 0)
	{
		printUsage(argv[0]);
		return -1;
	}

	if(cmdOptionExists(argv, argv+argc, "-p"))
	{
		parallel = true;
		char *nb_threads = getCmdOption(argv, argv+argc, "-p");

		if(nb_threads != NULL)
		{
			nb_parallel_threads = atoi(nb_threads);
			if(modulus_size_bits < 1)
			{
				printUsage(argv[0]);
				return -1;
			}
		}

		cout << "[[ Parallel ]]\t";
		cout << "NB THREADS ";
		if(nb_parallel_threads == 0)
			cout << "MAX_THREADS";
		else
			cout << nb_parallel_threads;
		cout << endl;
	}

	cout << endl << "[[ Using RSA modulus size "
		 <<	mpz_sizeinbase (N.get_mpz_t (), 2) << "]]" << endl;
	cout << "N = P*Q = " << N << endl << endl;

	//Initiate the Quadratic Sieve algorithm
	if(!parallel)
	{
		QS qs(N, NB_LINEAR_RELATIONS);

		//Start factoring
		qs.Factor();

		if(qs.GetFactor1 () == 1 || qs.GetFactor1 () == N)
		{
			cout << ">>>> Failed to factor " << N << " <<<<\t"
					<< "Try using more linear relations" << endl;
		}
		else
		{
			cout << endl << ">>>>>>> Factored " << N << endl;
			cout << "\t Factor 1: " << qs.GetFactor1 () << endl;
			cout << "\t Factor 2: " << qs.GetFactor2 () << endl;
		}
	}
	else
	{
		QSParallel qs(N, NB_LINEAR_RELATIONS);
		if(nb_parallel_threads != 0)
			qs.setNumThreads(nb_parallel_threads);

		//Start factoring
		qs.Factor();

		if(qs.GetFactor1 () == 1 || qs.GetFactor1 () == N)
		{
			cout << ">>>> Failed to factor " << N << " <<<<\t"
					<< "Try using more linear relations" << endl;
		}
		else
		{
			cout << endl << ">>>>>>> Factored " << N << endl;
			cout << "\t Factor 1: " << qs.GetFactor1 () << endl;
			cout << "\t Factor 2: " << qs.GetFactor2 () << endl;
		}
	}
	cout << endl << "Done!" << endl;
}
