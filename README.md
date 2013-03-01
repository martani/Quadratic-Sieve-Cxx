Quadratic Sieve C++
===================

A simple Quadratic Sieve implementation in C++, this is an adaptation of the C version here [https://github.com/martani/Quadratic-Sieve](https://github.com/martani/Quadratic-Sieve).

Status
------
The current version is fully working but have several performance issues compared to the `C` version.
It does not play nice with small RSA integers (less than `50 bits`) and natuarally with very large ones (more than `200 bits`).

Compilling and testing the code
-------------------
  1. `autoreconf -v`
  2. `./configure` (if another compiler, `clang` for example use: `CXX=clang ./configure`)
  3. `cd tests/ && make`
  4. `./test-qs`

Factoring RSA integers from C++ code
--------------------------------
You can use the code as follows:
  * Include the `QS.h` header file
  * Initialize the Qudratic Sieve algorithm from the `QS` class

````C++
    //N represents the RSA integer you want to factor
    //NB_LINEAR_RELATIONS is a const (between 1..10 usually) indicating how much it is likely to find a factorization
    QS qs(N, NB_LINEAR_RELATIONS);
````
  * Start factoring (this may take a while depending on the size of your RSA modulus)

````C++
    //Start factoring
    qs.Factor();
````
  * Get the results back

````C++
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
````

***Note*** If the algorithm fails to factor the RSA integer is returns the result `(N, 1)` as factors.

Check [tests/test-qs.C](tests/test-qs.C) for a complete example of how to use the `QS` class.