/*
 * matrix.h
 *
 *  Created on: Jan 17, 2013
 *      Author: sky
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include "gmp.h"
#include "../smooth-number.h"

/*
 * Represents a matrix using (gmp)mpz_t as rows
 */

class Matrix {
public:

	Matrix ()
	{
		this->n_rows = this->n_cols = 0;
		this->next_available_row = 0;
		_A = NULL;
	}

	Matrix (unsigned long int rows, unsigned long int cols)
	{
		if(rows < 1)
			return;

		this->n_rows = rows;
		this->n_cols = cols;

		this->next_available_row = 0;
		_A = new mpz_t [rows];
	}

	~Matrix ()
	{
		for(unsigned long int i=0; i<this->next_available_row; ++i)
			mpz_clear(this->_A[i]);

		delete [] this->_A;
	}

	void Init (unsigned long int rows, unsigned long int cols)
	{
		this->n_rows = rows;
		this->n_cols = cols;

		this->next_available_row = 0;

		//Release previous data is any
		if(this->_A != NULL)
		{
			for(unsigned long int i=0; i<this->next_available_row; ++i)
				mpz_clear(this->_A[i]);

			delete [] this->_A;
		}

		_A = new mpz_t [rows];
	}

	mpz_t &operator [] (unsigned long int row_index)
	{
		return _A[row_index];
	}

	const mpz_t &operator [] (unsigned long int row_index) const
	{
		return _A[row_index];
	}

	static Matrix *GetIdentity (unsigned long int nb_rows);

	//Push a new row at the end of the matrix.
	//If there are no more available rows, returns false, true otherwise.
	bool PushExponentVector (SmoothNumber& n);

	bool PushRow(const mpz_t r);

	unsigned long int row_dim () const { return this->n_rows; }
	unsigned long int col_dim () const { return this->n_cols; }



private:

	//The number of rows and columns in the matrix
	unsigned long int n_rows, n_cols;

	//The next free row in the matrix (goes from 0, 1, 2, ...)
	unsigned long int next_available_row;

	//The actual rows of type mpz_t
	mpz_t *_A;

	//The Identity matrix, used to keep track of linear relations between the rows
	//after the G

	Matrix(const Matrix& other);
	Matrix& operator= (const Matrix& other);
};


Matrix *Matrix::GetIdentity (unsigned long int nb_rows)
{
	Matrix *M = new Matrix (nb_rows, nb_rows);
	M->next_available_row = nb_rows;

	for(unsigned long int i=0; i<nb_rows; ++i)
	{
		//Init a nb_rows bit vector, set all entries to 0
		mpz_init2((*M)[i], nb_rows);

		//Set bit on diagonal at row i to 1
		mpz_setbit((*M)[i], i);
	}

	return M;
}


bool Matrix::PushExponentVector (SmoothNumber& n)
{
	//Is there more room to push the row?
	if(this->next_available_row >= this->n_rows)
		return false;

	mpz_init_set(this->_A[this->next_available_row], n.GetExponentVector ());
	this->next_available_row++;

	return true;
}

bool Matrix::PushRow(const mpz_t r)
{
	if(this->next_available_row >= this->n_rows)
		return false;

	mpz_init_set(this->_A[this->next_available_row], r);
	this->next_available_row++;

	return true;
}

#endif /* MATRIX_H_ */
