/*
 * gauss-elimination.h
 *
 *  Created on: Jan 17, 2013
 *      Author: sky
 */

#ifndef GAUSS_ELIMINATION_H_
#define GAUSS_ELIMINATION_H_

#include "gmp.h"
#include "matrix.h"

class GaussElimination {
private:

	//Fills the vector linear_relations after Echelonize is called
	void FillLinearRelations (Matrix& M, Matrix& ID);

	bool is_echelonized;

	mpz_t *linear_relations;
	int nb_linear_relations;

public:

	GaussElimination () : is_echelonized (false), linear_relations (NULL),
						nb_linear_relations (0) {}

	~GaussElimination ()
	{
		if(linear_relations != NULL)
			delete [] linear_relations;
	}

	//Performs a Gaussian elimination over the matrix M
	void Echelonize (Matrix& M);

	//Returns the linear relations resulted from the echelon form of the Matrix M
	const mpz_t* GetLinearRelations () const;
};


//Performs a Gaussian elimination over M
void GaussElimination::Echelonize (Matrix& M)
{
	Matrix *ID_p = Matrix::GetIdentity(M.row_dim ());
	Matrix &ID = *ID_p;

	unsigned long int current_non_reduced_row = 0;
	unsigned long int current_pivot_idx;
	unsigned long int next_pivot_idx;

	for(unsigned long int col=0; col<M.col_dim (); ++col)
	{
//		std::cout << "Column " << col << std::ends
//				<< " Current non reduced " << current_non_reduced_row
//				<< std::endl;

		next_pivot_idx = (unsigned long int)-1;

		//Get the first row which has the bit at column col equals to 1
		for(unsigned long int i=current_non_reduced_row; i<M.row_dim (); ++i)
		{
			//std::cout << "\ti " << i << std::endl;
			//std::cout << "M[i] = " << M[i] << std::endl;
			if(mpz_tstbit(M[i], col))
			{
				next_pivot_idx = i;
				break;
			}
		}

		//Now condidate row found, continue to next column
		if(next_pivot_idx == (unsigned long int)-1)
		{
//			std::cout << "No pivots at column " << col << std::endl;
			continue;
		}

		//The pivot is different from the current non reduced row, switch them!
		if(next_pivot_idx != current_non_reduced_row)
		{
			mpz_swap(M[next_pivot_idx], M[current_non_reduced_row]);
			mpz_swap(ID[next_pivot_idx], ID[current_non_reduced_row]);
		}
		current_pivot_idx = current_non_reduced_row;
		++current_non_reduced_row;

		//Reduce all the remaining rows with the pivot now located at
		for(unsigned long int i=current_non_reduced_row; i<M.row_dim (); ++i)
		{
			if(mpz_tstbit(M[i], col))
			{
				//XOR to eliminate the 1 in row i, column col.
				mpz_xor(M[i], M[i], M[current_pivot_idx]);
				mpz_xor(ID[i], ID[i], ID[current_pivot_idx]);
			}
		}
	}

	this->FillLinearRelations(M, ID);

	is_echelonized = true;
}


void GaussElimination::FillLinearRelations (Matrix& M, Matrix& ID)
{
	if(ID.row_dim() < 1)
		return;

	unsigned long int last_row_idx = ID.row_dim() - 1;
	int nb_linear_rel;

	while(mpz_cmp_ui(M[last_row_idx], 0) == 0)
	{
		--last_row_idx;
		++nb_linear_rel;
	}

	if(nb_linear_rel < 1)
		return;

	this->nb_linear_relations = nb_linear_rel;
	this->linear_relations = new mpz_t [this->nb_linear_relations];
	nb_linear_rel = 0;

	last_row_idx = ID.row_dim() - 1;
	while(mpz_cmp_ui(M[last_row_idx], 0) == 0)
	{
		mpz_init_set(this->linear_relations[nb_linear_rel], ID[last_row_idx]);
		--last_row_idx;
	}
}

const mpz_t* GaussElimination::GetLinearRelations () const
{
	if(!this->is_echelonized)
		throw std::logic_error ("Must Call Echelonize first");

	return this->linear_relations;
}

#endif /* GAUSS_ELIMINATION_H_ */
