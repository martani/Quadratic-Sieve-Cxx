/*
 * test-vector-smooth-number.C
 *
 *  Created on: Feb 28, 2013
 *      Author: sky
 */

#include <iostream>
#include <vector>
#include <omp.h>
#include "../src/smooth-number.h"

using namespace std;

/**
 * Check the result within valgrind to detect any memory leaks due to the
 * copy/assignment constructors.
 */

int main(int argc, char **argv)
{
	mpz_class x = 2;
	mpz_class N = 1234234234;

	vector<SmoothNumber> v;
	SmoothNumber tmp;
	v.push_back(tmp);
	v.back().InitWithoutExponentVector(x, N);

	SmoothNumber *a = new SmoothNumber (tmp);
	a->InitWithoutExponentVector(x, N);
	delete a;

	cout << "Done !" << endl;

	int nProcessors = omp_get_max_threads();
	cout << "OMP NB processors : " << nProcessors << endl;
}
