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

  //*a = static_cast<FLOAT *>(aligned_malloc(sizeof(FLOAT), ALIGN));
  *A = static_cast<FLOAT *>(aligned_malloc(2 * PARAM0*PARAM1*sizeof(FLOAT), ALIGN));
  *B = static_cast<FLOAT *>(aligned_malloc(2 * PARAM1*PARAM2*sizeof(FLOAT), ALIGN));
  //*b = static_cast<FLOAT *>(aligned_malloc(sizeof(FLOAT), ALIGN));
  *initC = static_cast<FLOAT *>(aligned_malloc(2 * PARAM0*PARAM2*sizeof(FLOAT), ALIGN));
  *C = static_cast<FLOAT *>(aligned_malloc(2 * PARAM0*PARAM2*sizeof(FLOAT), ALIGN));

  //rands(*a, 1, 1);
  rands(*A, PARAM0, 2 * PARAM1);
  rands(*B, PARAM1, 2 * PARAM2);
  //rands(*b, 1, 1);
  rands(*initC, PARAM0, 2 * PARAM2);
  memcpy(*C, *initC, 2 * PARAM0*PARAM2*sizeof(FLOAT));
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

  bool success = true;
  std::vector<FLOAT> temp(2*PARAM0*PARAM2, 0.);

  std::vector<string> errMsgs;

  for (int i = 0; i < PARAM0; ++i) {
      for (int j = 0; j < PARAM2; ++j) {
    	  for (int k = 0; k < PARAM1; ++k) {
              int creal_index = (2*PARAM2*i) + 8*(j/4) + (j % 4);
              int cimg_index = (2*PARAM2*i) + 8*(j/4) + (j % 4) + 4;

              int areal_index = (2*PARAM1*i) + 8*(k/4) + (k % 4);
              int aimg_index = (2*PARAM1*i) + 8*(k/4) + (k % 4) + 4;

              int breal_index = (2*PARAM2*k) + 8*(j/4) + (j % 4);
              int bimg_index = (2*PARAM2*k) + 8*(j/4) + (j % 4) + 4;
    		  temp[creal_index] += A[areal_index] * B[breal_index];
              temp[creal_index] -= A[aimg_index] * B[bimg_index];

              temp[cimg_index] += A[areal_index] * B[bimg_index];
              temp[cimg_index] += A[aimg_index] * B[breal_index];
    	  }
    	  //temp[i*PARAM2 + j] *= a;
      }
  }

  for (int i = 0; i < PARAM0; ++i) {
	  for (int j = 0; j < 2 * PARAM2; ++j) {
		  temp[2*i*PARAM2+j] += initC[2*i*PARAM2+j];

		  double err = fabs(C[2*i*PARAM2+j] - temp[2*i*PARAM2+j])/temp[2*i*PARAM2+j];
		  if(err > threshold)
		  {
			  success = false;
			  stringstream ss;
			  ss << "Error at (" << i << ", " << j << "): ";
			  ss << "C = " << C[2*i*PARAM2+j] << "\t-- Cref = " << temp[2*i*PARAM2+j] << "\t-- Err = " << err << endl;
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

  double flops = 4.*PARAM0*PARAM2*PARAM1;

  for (int k = 0; k < Rep; k++) {

	  memcpy(C, initC, 2 * PARAM0*PARAM2*sizeof(FLOAT));

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
