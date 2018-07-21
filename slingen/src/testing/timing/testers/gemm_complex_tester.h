/*
 * gemm_tester.h
 *
 *  Created on: June 6, 2012
 *      Author: danieles
 */

#pragma once

#include <iostream>
#include <sstream>
#include <ctime>
#include <list>
#include <cstdlib>
#include <cstring>
#include <cmath>

#ifdef TEST
#include "tsc.h"
#endif

#include "helpers.h"
#include "CommonDefs.h"

#include "kernels/gemm_complex_kernel.h"


/*
 * PARAM0 -> M
 * PARAM1 -> K
 * PARAM2 -> N
 */

inline void build(FLOAT ** A, FLOAT ** B,FLOAT ** initC, FLOAT ** C)
{
  srand(time(NULL));
  int col1 = PARAM1, col2 = PARAM2;
  if(PARAM1 % 4 != 0){
    col1 += (4 - (PARAM1 % 4));
  }
  if(PARAM2 % 4 != 0){
    col2 += (4 - (PARAM2 % 4));
  }
  //*a = static_cast<FLOAT *>(aligned_malloc(sizeof(FLOAT), ALIGN));
  *A = static_cast<FLOAT *>(aligned_malloc(2 * PARAM0*col1*sizeof(FLOAT), ALIGN));
  *B = static_cast<FLOAT *>(aligned_malloc(2 * PARAM1*col2*sizeof(FLOAT), ALIGN));
  //*b = static_cast<FLOAT *>(aligned_malloc(sizeof(FLOAT), ALIGN));
  *initC = static_cast<FLOAT *>(aligned_malloc(2 * PARAM0*col2*sizeof(FLOAT), ALIGN));
  *C = static_cast<FLOAT *>(aligned_malloc(2 * PARAM0*col2*sizeof(FLOAT), ALIGN));

  //rands(*a, 1, 1);
  comp_rands(*A, PARAM0, PARAM1);
  comp_rands(*B, PARAM1, PARAM2);
  //rands(*b, 1, 1);
  comp_rands(*initC, PARAM0, PARAM2);
  memcpy(*C, *initC, 2 * PARAM0*col2*sizeof(FLOAT));
}

inline void destroy(FLOAT * A, FLOAT * B,FLOAT * initC, FLOAT * C)
{
  //aligned_free(a);
  aligned_free(A);
  aligned_free(B);
  //aligned_free(b);
  aligned_free(initC);
  aligned_free(C);
}

int validate(FLOAT * A, FLOAT * B,FLOAT * initC, FLOAT * C, double threshold)
{
  int col1 = PARAM1, col2 = PARAM2;
  if(PARAM1 % 4 != 0){
    col1 += (4 - (PARAM1 % 4));
  }
  if(PARAM2 % 4 != 0){
    col2 += (4 - (PARAM2 % 4));
  }

  bool success = true;
  std::vector<FLOAT> temp(2*PARAM0*col2, 0.);

  std::vector<string> errMsgs;

  for (int i = 0; i < PARAM0; ++i) {
      for (int j = 0; j < PARAM2; ++j) {
    	  for (int k = 0; k < PARAM1; ++k) {
              int creal_index = (2*col2*i) + 8*(j/4) + (j % 4);
              int cimg_index = (2*col2*i) + 8*(j/4) + (j % 4) + 4;

              int areal_index = (2*col1*i) + 8*(k/4) + (k % 4);
              int aimg_index = (2*col1*i) + 8*(k/4) + (k % 4) + 4;

              int breal_index = (2*col2*k) + 8*(j/4) + (j % 4);
              int bimg_index = (2*col2*k) + 8*(j/4) + (j % 4) + 4;
             /* cout << "C_REAL INDEX : " << creal_index<<endl;
              cout << "C_IMG INDEX : " << cimg_index<<endl;
              cout << "A_REAL INDEX : " << areal_index<<endl;
              cout << "A_IMG INDEX : " << aimg_index<<endl;
              cout << "B_REAL INDEX : " << breal_index<<endl;
              cout << "B_IMG INDEX : " << bimg_index<<endl;
            */
    		  temp[creal_index] += A[areal_index] * B[breal_index];
              temp[creal_index] -= A[aimg_index] * B[bimg_index];

              temp[cimg_index] += A[areal_index] * B[bimg_index];
              temp[cimg_index] += A[aimg_index] * B[breal_index];
    	  }
    	  //temp[i*PARAM2 + j] *= a;
      }
  }

  for (int i = 0; i < PARAM0; ++i) {
	  for (int j = 0; j < PARAM2; ++j) {
          int creal_index = (2*col2*i) + 8*(j/4) + (j % 4);
		  temp[creal_index] += initC[creal_index];
		  //cout << "REAL INDEX : " << creal_index<<endl;

		  double err = fabs(C[creal_index] - temp[creal_index])/temp[creal_index];
		  if(err > threshold)
		  {
			  success = false;
			  stringstream ss;
			  ss << "Error at (" << i << ", " << j << ", real" << "): ";
			  ss << "C = " << C[creal_index] << "\t-- Cref = " << temp[creal_index] << "\t-- Err = " << err << endl;
			  errMsgs.push_back(ss.str());
		  }
	  }
  }
  for (int i = 0; i < PARAM0; ++i) {
	  for (int j = 0; j < PARAM2; ++j) {
          int cimg_index = (2*col2*i) + 8*(j/4) + (j % 4) + 4;
		  temp[cimg_index] += initC[cimg_index];
		  //cout << "IMG INDEX : " << cimg_index<<endl;

		  double err = fabs(C[cimg_index] - temp[cimg_index])/temp[cimg_index];
		  if(err > threshold)
		  {
			  success = false;
			  stringstream ss;
			  ss << "Error at (" << i << ", " << j << ", img" << "): ";
			  ss << "C = " << C[cimg_index] << "\t-- Cref = " << temp[cimg_index] << "\t-- Err = " << err << endl;
			  errMsgs.push_back(ss.str());
		  }
	  }
  }

  if(!success)
    for(std::vector<string>::const_iterator i = errMsgs.begin(); i != errMsgs.end(); i++)
      cout << *i;

  return !success;
}

/**
 * Test for gemm kernels
 */

int test()
{

  FLOAT * A, * B, * initC, * C;
  // myInt64 start, end, overhead;
  // double cycles = 0.;
  // size_t num_runs = RUNS, multiplier = 1;
  int retCode = 0;

  build(&A, &B, &initC, &C);
  int col1 = PARAM1, col2 = PARAM2;
  if(PARAM1 % 4 != 0){
    col1 += (4 - (PARAM1 % 4));
  }
  if(PARAM2 % 4 != 0){
    col2 += (4 - (PARAM2 % 4));
  }

#ifdef VALIDATE
  kernel(A,B,C);
  retCode = validate(A,B,initC,C, ERRTHRESH);
  if (retCode > 0) {
	cout << "Validation failed." << endl;
	return retCode;
  }
#endif

#ifdef TEST
  //Cache warm-up
  // RDTSCP reads ts register guaranteeing that the execution of all the code
  // we wanted to measure is completed. This way we avoid including the
  // execution of a CPUID in between. The last CPUID guarantees no other
  // instruction can be scheduled before it (and so also before RDTSCP)
  myInt64 start, end, overhead;
  double cycles = 0.;
  size_t num_runs = RUNS, multiplier = 1;

  init_tsc();
  overhead = get_tsc_overhead();

  do{
      num_runs = num_runs * multiplier;
      start = start_tsc();
      for(size_t i = 0; i < num_runs; i++) {
    	  kernel(A,B,C);
      }
      end = stop_tsc(start);
      if (end > overhead)
          end -= overhead;

      cycles = (double) end;
      multiplier = ceil (  (CYCLES_REQUIRED) / (cycles)  + 1.0 );

  }while (multiplier > 2);

  list< double > cycleList, flopList;
  size_t Rep = NUMREP;

  //the calculation for flops changes here
  double flops = 8.*PARAM0*PARAM2*PARAM1;

  for (int k = 0; k < Rep; k++) {

	  memcpy(C, initC, 2 * PARAM0*col2*sizeof(FLOAT));

      start = start_tsc();
      for (int i = 0; i < num_runs; ++i) {
    	  kernel(A,B,C);
      }
      end = stop_tsc(start);
      end -= overhead;

      cycles = ((double) end) / num_runs;

      cycleList.push_back(cycles);
      flopList.push_back(flops);

  }

  dumpList(cycleList, string(EXEC_PATH) + "/cycles.txt");
  dumpList(flopList, string(EXEC_PATH) + "/flops.txt");
#endif

  destroy(A, B,initC, C);

  return retCode;
}