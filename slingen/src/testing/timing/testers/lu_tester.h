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

#include "kernels/lu_kernel.h"


/*
 * PARAM0 -> M
 */

inline void build(FLOAT ** L, FLOAT ** U,FLOAT ** initC, FLOAT ** C)
{
  srand(time(NULL));

  //*a = static_cast<FLOAT *>(aligned_malloc(sizeof(FLOAT), ALIGN));
  *L = static_cast<FLOAT *>(aligned_malloc(PARAM0*PARAM0*sizeof(FLOAT), ALIGN));
  *U = static_cast<FLOAT *>(aligned_malloc(PARAM0*PARAM0*sizeof(FLOAT), ALIGN));
  //*b = static_cast<FLOAT *>(aligned_malloc(sizeof(FLOAT), ALIGN));
  *initC = static_cast<FLOAT *>(aligned_malloc(PARAM0*PARAM0*sizeof(FLOAT), ALIGN));
  *C = static_cast<FLOAT *>(aligned_malloc(PARAM0*PARAM0*sizeof(FLOAT), ALIGN));

  //rands(*a, 1, 1);
  lrands(*L, PARAM0, PARAM0);
  urands(*U, PARAM0, PARAM0);
  //rands(*b, 1, 1);
  rands(*initC, PARAM0, PARAM0);
  memcpy(*C, *initC, PARAM0*PARAM0*sizeof(FLOAT));
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
  std::vector<FLOAT> temp(PARAM0*PARAM0, 0.);

  std::vector<string> errMsgs;

  for (int i = 0; i < PARAM0; ++i) {
      for (int j = 0; j < PARAM0; ++j) {
    	  for (int k = 0; k < PARAM0; ++k) {
    		  temp[i*PARAM0 + j] += A[i*PARAM0 + k] * B[k*PARAM0 + j];
    	  }
    	  //temp[i*PARAM2 + j] *= a;
      }
  }

  for (int i = 0; i < PARAM0; ++i) {
	  for (int j = 0; j < PARAM0; ++j) {
		  temp[i*PARAM0+j] += initC[i*PARAM0+j];

		  double err = fabs(C[i*PARAM0+j] - temp[i*PARAM0+j])/temp[i*PARAM0+j];
		  if(err > threshold)
		  {
			  success = false;
			  stringstream ss;
			  ss << "Error at (" << i << ", " << j << "): ";
			  ss << "C = " << C[i*PARAM0+j] << "\t-- Cref = " << temp[i*PARAM0+j] << "\t-- Err = " << err << endl;
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

  double flops = 2.*PARAM0*PARAM0*PARAM0;

  for (int k = 0; k < Rep; k++) {

	  memcpy(C, initC, PARAM0*PARAM0*sizeof(FLOAT));

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
