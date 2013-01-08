/*
 * utils.h
 *
 *  Created on: Jan 5, 2013
 *      Author: martani
 *      Copyright martani 2013
 */

#include <sys/time.h>

#ifndef UTILS_H_
#define UTILS_H_


#define TIMER_DECLARE(s)  		\
	timeval __timer_q_##s; 		\
	timeval __timer_q_end_##s;	\
	double  __timer_q_diff_##s = 0.0;

#define TIMER_START(s)    gettimeofday(&__timer_q_##s, NULL);

#define TIMER_STOP(s) 	  							\
	gettimeofday(&__timer_q_end_##s, NULL);			\
	__timer_q_diff_##s += (__timer_q_end_##s.tv_sec  - __timer_q_##s.tv_sec) * 1000.0;	\
	__timer_q_diff_##s += (__timer_q_end_##s.tv_usec - __timer_q_##s.tv_usec) / 1000.0;

#define TIMER_REPORT(s)													\
	std::cout << "Total [" #s "] time: " << __timer_q_diff_##s << " ms."	\
			  << std::endl;


#endif /* UTILS_H_ */
