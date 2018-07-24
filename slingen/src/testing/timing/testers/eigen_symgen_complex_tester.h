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
#include<complex>

#ifdef TEST
#include "tsc.h"
#endif

#include<Eigen/Dense>
using namespace Eigen;
using namespace std;
#include "helpers.h"
#include"params.h"
#include "CommonDefs.h"

#define MATRIX vector<complex<FLOAT> >

#define MAT Map< MatrixXcd, RowMajor >

/*
 * PARAM0 -> M
 * PARAM1 -> K
 * PARAM2 -> N
 */
/*
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
  rands(*A, PARAM0, 2*PARAM1);
  rands(*B, PARAM1, 2*PARAM2);
  //rands(*b, 1, 1);
  rands(*initC, PARAM0, 2*PARAM2);
  memcpy(*C, *initC, 2 * PARAM0*PARAM2*sizeof(FLOAT));
}
/*
inline void destroy(FLOAT * A, FLOAT * B,FLOAT * initC, FLOAT * C)
{
  //aligned_free(a);
  aligned_free(A);
  aligned_free(B);
  //aligned_free(b);
  aligned_free(initC);
  aligned_free(C);
}
*/
static __attribute__((noinline)) void kernel(MAT &A, MAT &B, MAT &C){
    C = (A*B);
}

int validate(MATRIX A, MATRIX B, MATRIX C, double threshold)
{

  bool success = true;
  std::vector<complex<FLOAT> > temp(PARAM0*PARAM0, (0.,0.));

  std::vector<string> errMsgs;

  for (int i = 0; i < PARAM0; ++i) {
         for (int j = 0; j < PARAM0; ++j) {
             int cindex = PARAM0*j + i;
             for (int k = 0; k < PARAM0; ++k) {

                 int aindex = PARAM0*k + i;

                 int bindex = PARAM0*j + k;

                 temp[cindex] += A[aindex] * B[bindex];

             }
       //          cout<<C[cindex]<<endl;
             //temp[i*PARAM2 + j] *= a;
         }
     }

  for (int i = 0; i < PARAM0; ++i) {
	  for (int j = 0; j < PARAM0; ++j) {
		  //temp[i*PARAM0+j].real() += initC[i*PARAM0+j].real();

		  double err = fabs(C[i*PARAM0+j].real() - temp[i*PARAM0+j].real())/temp[i*PARAM0+j].real();
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
  for (int i = 0; i < PARAM0; ++i) {
	  for (int j = 0; j < PARAM0; ++j) {
		  //temp[i*PARAM2+j].imag() += initC[i*PARAM2+j].imag();

		  double err = fabs(C[i*PARAM0+j].imag() - temp[i*PARAM0+j].imag())/temp[i*PARAM0+j].imag();
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

  MATRIX A, B, C;
  // myInt64 start, end, overhead;
  // double cycles = 0.;
  // size_t num_runs = RUNS, multiplier = 1;
  int retCode = 0;
  A.resize(PARAM0*PARAM0);
  B.resize(PARAM0*PARAM0);
  C.resize(PARAM0*PARAM0);
  //build(&A, &B, &initC, &C);
  for(int i=0;i<PARAM0;++i){
    for(int j=i;j<PARAM0;++j){
    complex<FLOAT> temp((FLOAT)(rand())/RAND_MAX + 1 ,(FLOAT)(rand())/RAND_MAX + 1);
    A[j*PARAM0 + i] = A[i*PARAM0 + j] = temp;

    }
  }
  for(int i=0;i<PARAM0*PARAM0;++i){
    complex<FLOAT> temp((FLOAT)(rand())/RAND_MAX + 1,(FLOAT)(rand())/RAND_MAX + 1);
    B[i] = temp;
  }
  MAT aMat(&A[0] , PARAM0 , PARAM0);
  MAT bMat(&B[0] , PARAM0 , PARAM0);
  MAT cMat(&C[0] , PARAM0 , PARAM0);
#ifdef VALIDATE
  kernel(aMat,bMat,cMat);
  retCode = validate(A,B,C, ERRTHRESH);
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
    	  kernel(aMat,bMat,cMat);
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
  double flops = 6.*PARAM0*PARAM0*PARAM0;

  for (int k = 0; k < Rep; k++) {

	  //memcpy(C, initC, 2 * PARAM0*PARAM2*sizeof(FLOAT));
      //C = initC;
      //MAT cMat(&C[0],PARAM0,PARAM2);
      start = start_tsc();
      for (int i = 0; i < num_runs; ++i) {
    	  kernel(aMat,bMat,cMat);
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

  //destroy(A, B,initC, C);

  return retCode;
}