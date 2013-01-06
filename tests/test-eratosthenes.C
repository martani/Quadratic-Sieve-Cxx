#include <iostream>
#include "../src/eratosthenes.h"
#include "../src/utils.h"

using namespace std;

int main(void)
{

	uint64_t base = 1000*1000*100;
	vector<uint64_t> primes;
	Erastosthenes eratosthenesSieve (false);

	TIMER_DECLARE(timer)
	TIMER_START(timer)
	eratosthenesSieve.GetPrimes(primes, base);
	TIMER_STOP(timer)

	cout << "Found " << primes.size() << " primes (<= " << base << ")" << endl;
	TIMER_REPORT(timer)

	mpz_class N;
	N.set_str("45342314233", 10);

	TIMER_DECLARE(timer2)
	TIMER_START(timer2)
	eratosthenesSieve.GetPrimes_QuadraticResidue(primes, base, N);
	TIMER_STOP(timer2)

	cout << "Found " << primes.size() << " primes (<= " << base << ")" << endl;
	TIMER_REPORT(timer2)

	return 0;
}
