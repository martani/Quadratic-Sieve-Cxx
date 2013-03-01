/*
 * utils.h
 *
 *  Created on: Jan 5, 2013
 *      Author: martani
 *      Copyright martani 2013
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <sys/time.h>
#include <iostream>
#include <vector>
#include <fstream>

#include "debug-cout.h"
#include "linear-algebra/matrix.h"


#ifndef	ENABLE_TIMERS
#define TIMER_DECLARE(s)  		\
	timeval __timer_q_##s; 		\
	timeval __timer_q_end_##s;	\
	double  __timer_q_diff_##s = 0.0;

#define TIMER_START(s)    gettimeofday(&__timer_q_##s, NULL);

#define TIMER_STOP(s) 	  							\
	gettimeofday(&__timer_q_end_##s, NULL);			\
	__timer_q_diff_##s += (__timer_q_end_##s.tv_sec  - __timer_q_##s.tv_sec) * 1000.0;	\
	__timer_q_diff_##s += (__timer_q_end_##s.tv_usec - __timer_q_##s.tv_usec) / 1000.0;

#define TIMER_REPORT(s)			\
	std::cout << "Total [" #s "] time: " << __timer_q_diff_##s << " ms."	\
			  << std::endl;
#else
#define TIMER_DECLARE(s)
#define TIMER_START(s)
#define TIMER_STOP(s)
#define TIMER_REPORT(s)
#endif

namespace Utils {

/*
 * Appends a tab \t separated values of the vector v to a file.
 * Inserts a new line (\t\n) at the end.
 */
template<typename T>
static void AppendVectorToFile(vector<T>& v, std::string file_name)
{
	ofstream file;
	file.open (file_name.c_str (), ios::app);

	for(int i=0; i<v.size (); ++i)
	{
		file << v[i] << "\t";
	}
	file << "\n";

	file.close ();
}



static void dumpMatrixAsPbmImage(const Matrix& A, std::string output_file_name)
{
	char buffer[512];
	unsigned char output_byte = 0;

	FILE *outStream = fopen(output_file_name.c_str (), "wb");

	//magic PBM header
#ifdef __LP64__	//64 bit machine
	sprintf(buffer, "P4\n# matrix size(%lu, %lu)\n%lu %lu\n", A.row_dim (), A.col_dim (), A.col_dim (), A.row_dim ());
#else			//32 bit machine
	sprintf(buffer, "P4\n# matrix size(%u, %u)\n%u %u\n", A.row_dim (), A.col_dim (), A.coldim (), A.rowdim ());
#endif

	fwrite(buffer, sizeof(char), strlen(buffer), outStream);


	for(int row=0; row<A.row_dim(); ++row)
	{
		int col;
		//for each bit in the mpz_t row
		for(col=0; col<A.col_dim(); ++col)
		{
			if(mpz_tstbit(A[row], col) == 1)
			{
				output_byte |= (1 << (7 - (col%8)));
			}
			else
			{
				output_byte &= ~(1 << (7 - (col%8)));
			}

			if(col%8 == 7) //flush byte every 8 cols
			{
				fwrite(&output_byte, sizeof(unsigned char), 1, outStream);
				output_byte = 0;
			}
		}

		if(col%8 != 0)
			fwrite(&output_byte, sizeof(unsigned char), 1, outStream);

		fflush(outStream);
	}

	fclose(outStream);
//		row_size = i_A->size ();
//
//		j=0;
//		for(col=0; col<m; ++col){
//
//			if(j<row_size && (*i_A)[j].first == col)
//			{
//				output_byte |= (1 << (7 - (col%8)));
//				j++;
//			}
//			else
//			{
//				output_byte &= ~(1 << (7 - (col%8)));
//			}
//
//			if(col%8 == 7) //flush byte every 8 cols
//			{
//				fwrite(&output_byte, sizeof(unsigned char), 1, outStream);
//				output_byte = 0;
//			}
//		}
//
//		if(col%8 != 0)
//			fwrite(&output_byte, sizeof(unsigned char), 1, outStream);
//
//		fflush(outStream);
//
//		++i_A;
//	}
//
//	fclose(outStream);
}


}

#endif /* UTILS_H_ */
