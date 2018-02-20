/*
 * ltrtri_kernel.h
 *
Decl { {u'L': LowerTriangular[L, (28, 28), GenMatAccess]} }
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ann: {'part_schemes': {'tril_inv_ow_opt': {'m': 'm1.ll'}}, 'cl1ck_v': 1, 'variant_tag': 'tril_inv_ow_opt_m1'}

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Entry 0:
Eq: Tile( (1, 1), G(h(1, 28, 0), L[28,28],h(1, 28, 0)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, 0), L[28,28],h(1, 28, 0)) ) )
Eq.ann: {}
Entry 1:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(3, 28, 1), L[28,28],h(1, 28, 0)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(3, 28, 1), L[28,28],h(1, 28, 0)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 0), L[28,28],h(1, 28, 0)) ) ) )
Eq.ann: {}
Entry 2:
Eq: Tile( (1, 1), G(h(1, 28, 1), L[28,28],h(1, 28, 1)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, 1), L[28,28],h(1, 28, 1)) ) )
Eq.ann: {}
Entry 3:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 2), L[28,28],h(1, 28, 1)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 2), L[28,28],h(1, 28, 1)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 1), L[28,28],h(1, 28, 1)) ) ) )
Eq.ann: {}
Entry 4:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 2), L[28,28],h(1, 28, 0)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 2), L[28,28],h(1, 28, 0)) ) ) - ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 2), L[28,28],h(1, 28, 1)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 1), L[28,28],h(1, 28, 0)) ) ) ) )
Eq.ann: {}
Entry 5:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 1), L[28,28],h(1, 28, 0)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 1), L[28,28],h(1, 28, 1)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 1), L[28,28],h(1, 28, 0)) ) ) ) )
Eq.ann: {}
Entry 6:
Eq: Tile( (1, 1), G(h(1, 28, 2), L[28,28],h(1, 28, 2)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, 2), L[28,28],h(1, 28, 2)) ) )
Eq.ann: {}
Entry 7:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 3), L[28,28],h(1, 28, 2)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 3), L[28,28],h(1, 28, 2)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 2), L[28,28],h(1, 28, 2)) ) ) )
Eq.ann: {}
Entry 8:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 3), L[28,28],h(2, 28, 0)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 3), L[28,28],h(2, 28, 0)) ) ) - ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 3), L[28,28],h(1, 28, 2)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 2), L[28,28],h(2, 28, 0)) ) ) ) )
Eq.ann: {}
Entry 9:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 2), L[28,28],h(2, 28, 0)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 2), L[28,28],h(1, 28, 2)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 2), L[28,28],h(2, 28, 0)) ) ) ) )
Eq.ann: {}
Entry 10:
Eq: Tile( (1, 1), G(h(1, 28, 3), L[28,28],h(1, 28, 3)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, 3), L[28,28],h(1, 28, 3)) ) )
Eq.ann: {}
Entry 11:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 3), L[28,28],h(3, 28, 0)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 3), L[28,28],h(1, 28, 3)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 3), L[28,28],h(3, 28, 0)) ) ) ) )
Eq.ann: {}
Entry 12:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(24, 28, 4), L[28,28],h(4, 28, 0)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(24, 28, 4), L[28,28],h(4, 28, 0)) ) ) * Tile( (1, 1), Tile( (4, 4), G(h(4, 28, 0), L[28,28],h(4, 28, 0)) ) ) )
Eq.ann: {}
Entry 13:
For_{fi210;4;23;4} ( Entry 0:
Eq: Tile( (1, 1), G(h(1, 28, fi210), L[28,28],h(1, 28, fi210)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, fi210), L[28,28],h(1, 28, fi210)) ) )
Eq.ann: {}
Entry 1:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(3, 28, fi210 + 1), L[28,28],h(1, 28, fi210)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(3, 28, fi210 + 1), L[28,28],h(1, 28, fi210)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210), L[28,28],h(1, 28, fi210)) ) ) )
Eq.ann: {}
Entry 2:
Eq: Tile( (1, 1), G(h(1, 28, fi210 + 1), L[28,28],h(1, 28, fi210 + 1)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, fi210 + 1), L[28,28],h(1, 28, fi210 + 1)) ) )
Eq.ann: {}
Entry 3:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(2, 28, fi210 + 2), L[28,28],h(1, 28, fi210 + 1)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, fi210 + 2), L[28,28],h(1, 28, fi210 + 1)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 1), L[28,28],h(1, 28, fi210 + 1)) ) ) )
Eq.ann: {}
Entry 4:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(2, 28, fi210 + 2), L[28,28],h(1, 28, fi210)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, fi210 + 2), L[28,28],h(1, 28, fi210)) ) ) - ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, fi210 + 2), L[28,28],h(1, 28, fi210 + 1)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 1), L[28,28],h(1, 28, fi210)) ) ) ) )
Eq.ann: {}
Entry 5:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 1), L[28,28],h(1, 28, fi210)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 1), L[28,28],h(1, 28, fi210 + 1)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 1), L[28,28],h(1, 28, fi210)) ) ) ) )
Eq.ann: {}
Entry 6:
Eq: Tile( (1, 1), G(h(1, 28, fi210 + 2), L[28,28],h(1, 28, fi210 + 2)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, fi210 + 2), L[28,28],h(1, 28, fi210 + 2)) ) )
Eq.ann: {}
Entry 7:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 3), L[28,28],h(1, 28, fi210 + 2)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 3), L[28,28],h(1, 28, fi210 + 2)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 2), L[28,28],h(1, 28, fi210 + 2)) ) ) )
Eq.ann: {}
Entry 8:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 3), L[28,28],h(2, 28, fi210)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 3), L[28,28],h(2, 28, fi210)) ) ) - ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 3), L[28,28],h(1, 28, fi210 + 2)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 2), L[28,28],h(2, 28, fi210)) ) ) ) )
Eq.ann: {}
Entry 9:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 2), L[28,28],h(2, 28, fi210)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 2), L[28,28],h(1, 28, fi210 + 2)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 2), L[28,28],h(2, 28, fi210)) ) ) ) )
Eq.ann: {}
Entry 10:
Eq: Tile( (1, 1), G(h(1, 28, fi210 + 3), L[28,28],h(1, 28, fi210 + 3)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, fi210 + 3), L[28,28],h(1, 28, fi210 + 3)) ) )
Eq.ann: {}
Entry 11:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 3), L[28,28],h(3, 28, fi210)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 3), L[28,28],h(1, 28, fi210 + 3)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, fi210 + 3), L[28,28],h(3, 28, fi210)) ) ) ) )
Eq.ann: {}
Entry 12:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(-fi210 + 24, 28, fi210 + 4), L[28,28],h(4, 28, fi210)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(-fi210 + 24, 28, fi210 + 4), L[28,28],h(4, 28, fi210)) ) ) * Tile( (1, 1), Tile( (4, 4), G(h(4, 28, fi210), L[28,28],h(4, 28, fi210)) ) ) )
Eq.ann: {}
Entry 13:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(-fi210 + 24, 28, fi210 + 4), L[28,28],h(fi210, 28, 0)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(-fi210 + 24, 28, fi210 + 4), L[28,28],h(fi210, 28, 0)) ) ) - ( Tile( (1, 1), Tile( (4, 4), G(h(-fi210 + 24, 28, fi210 + 4), L[28,28],h(4, 28, fi210)) ) ) * Tile( (1, 1), Tile( (4, 4), G(h(4, 28, fi210), L[28,28],h(fi210, 28, 0)) ) ) ) )
Eq.ann: {}
Entry 14:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(4, 28, fi210), L[28,28],h(fi210, 28, 0)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(4, 28, fi210), L[28,28],h(4, 28, fi210)) ) ) * Tile( (1, 1), Tile( (4, 4), G(h(4, 28, fi210), L[28,28],h(fi210, 28, 0)) ) ) ) )
Eq.ann: {}
 )Entry 14:
Eq: Tile( (1, 1), G(h(1, 28, 24), L[28,28],h(1, 28, 24)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, 24), L[28,28],h(1, 28, 24)) ) )
Eq.ann: {}
Entry 15:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(3, 28, 25), L[28,28],h(1, 28, 24)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(3, 28, 25), L[28,28],h(1, 28, 24)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 24), L[28,28],h(1, 28, 24)) ) ) )
Eq.ann: {}
Entry 16:
Eq: Tile( (1, 1), G(h(1, 28, 25), L[28,28],h(1, 28, 25)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, 25), L[28,28],h(1, 28, 25)) ) )
Eq.ann: {}
Entry 17:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 26), L[28,28],h(1, 28, 25)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 26), L[28,28],h(1, 28, 25)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 25), L[28,28],h(1, 28, 25)) ) ) )
Eq.ann: {}
Entry 18:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 26), L[28,28],h(1, 28, 24)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 26), L[28,28],h(1, 28, 24)) ) ) - ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 26), L[28,28],h(1, 28, 25)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 25), L[28,28],h(1, 28, 24)) ) ) ) )
Eq.ann: {}
Entry 19:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 25), L[28,28],h(1, 28, 24)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 25), L[28,28],h(1, 28, 25)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 25), L[28,28],h(1, 28, 24)) ) ) ) )
Eq.ann: {}
Entry 20:
Eq: Tile( (1, 1), G(h(1, 28, 26), L[28,28],h(1, 28, 26)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, 26), L[28,28],h(1, 28, 26)) ) )
Eq.ann: {}
Entry 21:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 27), L[28,28],h(1, 28, 26)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 27), L[28,28],h(1, 28, 26)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 26), L[28,28],h(1, 28, 26)) ) ) )
Eq.ann: {}
Entry 22:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 27), L[28,28],h(2, 28, 24)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 27), L[28,28],h(2, 28, 24)) ) ) - ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 27), L[28,28],h(1, 28, 26)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 26), L[28,28],h(2, 28, 24)) ) ) ) )
Eq.ann: {}
Entry 23:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 26), L[28,28],h(2, 28, 24)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 26), L[28,28],h(1, 28, 26)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 26), L[28,28],h(2, 28, 24)) ) ) ) )
Eq.ann: {}
Entry 24:
Eq: Tile( (1, 1), G(h(1, 28, 27), L[28,28],h(1, 28, 27)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, 27), L[28,28],h(1, 28, 27)) ) )
Eq.ann: {}
Entry 25:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 27), L[28,28],h(3, 28, 24)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 27), L[28,28],h(1, 28, 27)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 27), L[28,28],h(3, 28, 24)) ) ) ) )
Eq.ann: {}
Entry 26:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(4, 28, 24), L[28,28],h(24, 28, 0)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(4, 28, 24), L[28,28],h(4, 28, 24)) ) ) * Tile( (1, 1), Tile( (4, 4), G(h(4, 28, 24), L[28,28],h(24, 28, 0)) ) ) ) )
Eq.ann: {}
 *
 * Created on: 2016-10-12
 * Author: danieles
 */

#pragma once

#include <x86intrin.h>


#define PARAM0 28

#define ERRTHRESH 1e-14

#define NUMREP 30

#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))


static __attribute__((noinline)) void kernel(double * L)
{
  __m256d _t0_0, _t0_1, _t0_2, _t0_3, _t0_4, _t0_5, _t0_6, _t0_7,
	_t0_8, _t0_9, _t0_10, _t0_11, _t0_12, _t0_13, _t0_14, _t0_15,
	_t0_16, _t0_17, _t0_18, _t0_19, _t0_20, _t0_21, _t0_22, _t0_23,
	_t0_24, _t0_25, _t0_26, _t0_27, _t0_28, _t0_29, _t0_30, _t0_31,
	_t0_32, _t0_33, _t0_34, _t0_35, _t0_36, _t0_37, _t0_38, _t0_39,
	_t0_40, _t0_41, _t0_42, _t0_43, _t0_44, _t0_45, _t0_46, _t0_47,
	_t0_48, _t0_49, _t0_50, _t0_51, _t0_52, _t0_53, _t0_54, _t0_55,
	_t0_56, _t0_57, _t0_58;
  __m256d _t1_0, _t1_1, _t1_2, _t1_3, _t1_4, _t1_5, _t1_6, _t1_7,
	_t1_8, _t1_9, _t1_10, _t1_11, _t1_12, _t1_13, _t1_14, _t1_15,
	_t1_16, _t1_17, _t1_18, _t1_19;
  __m256d _t2_0, _t2_1, _t2_2, _t2_3, _t2_4, _t2_5, _t2_6, _t2_7,
	_t2_8, _t2_9, _t2_10, _t2_11, _t2_12, _t2_13, _t2_14, _t2_15,
	_t2_16, _t2_17, _t2_18, _t2_19, _t2_20, _t2_21, _t2_22, _t2_23,
	_t2_24, _t2_25, _t2_26, _t2_27, _t2_28, _t2_29, _t2_30, _t2_31,
	_t2_32, _t2_33, _t2_34, _t2_35, _t2_36, _t2_37, _t2_38, _t2_39,
	_t2_40, _t2_41, _t2_42, _t2_43, _t2_44, _t2_45, _t2_46, _t2_47,
	_t2_48, _t2_49, _t2_50, _t2_51, _t2_52, _t2_53, _t2_54, _t2_55,
	_t2_56, _t2_57, _t2_58;
  __m256d _t3_0, _t3_1, _t3_2, _t3_3, _t3_4, _t3_5, _t3_6, _t3_7,
	_t3_8, _t3_9, _t3_10, _t3_11, _t3_12, _t3_13, _t3_14, _t3_15,
	_t3_16, _t3_17, _t3_18, _t3_19, _t3_20, _t3_21, _t3_22, _t3_23;
  __m256d _t4_0, _t4_1, _t4_2, _t4_3, _t4_4, _t4_5, _t4_6, _t4_7,
	_t4_8, _t4_9, _t4_10, _t4_11, _t4_12, _t4_13, _t4_14, _t4_15,
	_t4_16, _t4_17, _t4_18, _t4_19, _t4_20, _t4_21, _t4_22, _t4_23,
	_t4_24, _t4_25, _t4_26, _t4_27;
  __m256d _t5_0, _t5_1, _t5_2, _t5_3;
  __m256d _t6_0, _t6_1, _t6_2, _t6_3, _t6_4, _t6_5, _t6_6, _t6_7,
	_t6_8, _t6_9, _t6_10, _t6_11;
  __m256d _t7_0, _t7_1, _t7_2, _t7_3, _t7_4, _t7_5, _t7_6, _t7_7,
	_t7_8, _t7_9, _t7_10, _t7_11, _t7_12, _t7_13, _t7_14, _t7_15,
	_t7_16, _t7_17, _t7_18, _t7_19, _t7_20, _t7_21, _t7_22, _t7_23,
	_t7_24, _t7_25, _t7_26, _t7_27, _t7_28, _t7_29, _t7_30, _t7_31,
	_t7_32, _t7_33, _t7_34, _t7_35, _t7_36, _t7_37, _t7_38, _t7_39,
	_t7_40, _t7_41, _t7_42, _t7_43, _t7_44, _t7_45, _t7_46, _t7_47,
	_t7_48, _t7_49, _t7_50, _t7_51, _t7_52, _t7_53, _t7_54, _t7_55,
	_t7_56, _t7_57, _t7_58;
  __m256d _t8_0, _t8_1, _t8_2, _t8_3, _t8_4, _t8_5, _t8_6, _t8_7;

  _t0_0 = _mm256_maskload_pd(L, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
  _t0_1 = _mm256_permute2f128_pd(_mm256_unpacklo_pd(_mm256_maskload_pd(L + 28, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0)), _mm256_maskload_pd(L + 56, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0))), _mm256_maskload_pd(L + 84, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0)), 32);
  _t0_2 = _mm256_maskload_pd(L + 29, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
  _t0_3 = _mm256_shuffle_pd(_mm256_maskload_pd(L + 57, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0)), _mm256_maskload_pd(L + 85, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0)), 0);
  _t0_6 = _mm256_maskload_pd(L + 58, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
  _t0_7 = _mm256_maskload_pd(L + 86, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
  _t0_10 = _mm256_maskload_pd(L + 87, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));

  // Constant 1x1 -> 1x4
  _t0_44 = _mm256_set_pd(0, 0, 0, 1);

  // 1x1 -> 1x4
  _t0_49 = _t0_0;

  // 4-BLAC: 1x4 / 1x4
  _t0_18 = _mm256_div_pd(_t0_44, _t0_49);
  _t0_0 = _t0_18;

  // 3x1 -> 4x1
  _t0_19 = _t0_1;

  // 1x1 -> 1x4
  _t0_20 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_0, _t0_0, 32), _mm256_permute2f128_pd(_t0_0, _t0_0, 32), 0);

  // 4-BLAC: 4x1 Kro 1x4
  _t0_21 = _mm256_mul_pd(_t0_19, _t0_20);
  _t0_1 = _t0_21;

  // Constant 1x1 -> 1x4
  _t0_22 = _mm256_set_pd(0, 0, 0, 1);

  // 1x1 -> 1x4
  _t0_23 = _t0_2;

  // 4-BLAC: 1x4 / 1x4
  _t0_24 = _mm256_div_pd(_t0_22, _t0_23);
  _t0_2 = _t0_24;

  // 2x1 -> 4x1
  _t0_25 = _t0_3;

  // 1x1 -> 1x4
  _t0_26 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_2, _t0_2, 32), _mm256_permute2f128_pd(_t0_2, _t0_2, 32), 0);

  // 4-BLAC: 4x1 Kro 1x4
  _t0_27 = _mm256_mul_pd(_t0_25, _t0_26);
  _t0_3 = _t0_27;

  // 2x1 -> 4x1
  _t0_28 = _mm256_shuffle_pd(_mm256_blend_pd(_mm256_setzero_pd(), _t0_1, 2), _mm256_permute2f128_pd(_t0_1, _t0_1, 129), 5);

  // 2x1 -> 4x1
  _t0_29 = _t0_3;

  // 1x1 -> 1x4
  _t0_30 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_1, _t0_1, 32), _mm256_permute2f128_pd(_t0_1, _t0_1, 32), 0);

  // 4-BLAC: 4x1 Kro 1x4
  _t0_31 = _mm256_mul_pd(_t0_29, _t0_30);

  // 4-BLAC: 4x1 - 4x1
  _t0_32 = _mm256_sub_pd(_t0_28, _t0_31);
  _t0_4 = _t0_32;

  // 1x1 -> 1x4
  _t0_33 = _t0_2;

  // 1x1 -> 1x4
  _t0_34 = _mm256_blend_pd(_mm256_setzero_pd(), _t0_1, 1);

  // 4-BLAC: 1x4 Kro 1x4
  _t0_35 = _mm256_mul_pd(_t0_33, _t0_34);

  // 4-BLAC: -( 1x4 )
  _t0_36 = _mm256_sub_pd(_mm256_setzero_pd(), _t0_35);
  _t0_5 = _t0_36;

  // Constant 1x1 -> 1x4
  _t0_37 = _mm256_set_pd(0, 0, 0, 1);

  // 1x1 -> 1x4
  _t0_38 = _t0_6;

  // 4-BLAC: 1x4 / 1x4
  _t0_39 = _mm256_div_pd(_t0_37, _t0_38);
  _t0_6 = _t0_39;

  // 1x1 -> 1x4
  _t0_40 = _t0_7;

  // 1x1 -> 1x4
  _t0_41 = _t0_6;

  // 4-BLAC: 1x4 Kro 1x4
  _t0_42 = _mm256_mul_pd(_t0_40, _t0_41);
  _t0_7 = _t0_42;

  // 1x2 -> 1x4
  _t0_43 = _mm256_unpackhi_pd(_mm256_blend_pd(_t0_4, _mm256_setzero_pd(), 12), _mm256_blend_pd(_t0_3, _mm256_setzero_pd(), 12));

  // 1x1 -> 1x4
  _t0_45 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_7, _t0_7, 32), _mm256_permute2f128_pd(_t0_7, _t0_7, 32), 0);

  // 1x2 -> 1x4
  _t0_46 = _mm256_blend_pd(_mm256_unpacklo_pd(_t0_4, _t0_3), _mm256_setzero_pd(), 12);

  // 4-BLAC: 1x4 Kro 1x4
  _t0_47 = _mm256_mul_pd(_t0_45, _t0_46);

  // 4-BLAC: 1x4 - 1x4
  _t0_48 = _mm256_sub_pd(_t0_43, _t0_47);
  _t0_8 = _t0_48;

  // 1x1 -> 1x4
  _t0_50 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_6, _t0_6, 32), _mm256_permute2f128_pd(_t0_6, _t0_6, 32), 0);

  // 1x2 -> 1x4
  _t0_51 = _mm256_blend_pd(_mm256_unpacklo_pd(_t0_4, _t0_3), _mm256_setzero_pd(), 12);

  // 4-BLAC: 1x4 Kro 1x4
  _t0_52 = _mm256_mul_pd(_t0_50, _t0_51);

  // 4-BLAC: -( 1x4 )
  _t0_53 = _mm256_sub_pd(_mm256_setzero_pd(), _t0_52);
  _t0_9 = _t0_53;

  // Constant 1x1 -> 1x4
  _t0_54 = _mm256_set_pd(0, 0, 0, 1);

  // 1x1 -> 1x4
  _t0_55 = _t0_10;

  // 4-BLAC: 1x4 / 1x4
  _t0_56 = _mm256_div_pd(_t0_54, _t0_55);
  _t0_10 = _t0_56;

  // 1x1 -> 1x4
  _t0_57 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_10, _t0_10, 32), _mm256_permute2f128_pd(_t0_10, _t0_10, 32), 0);

  // 1x3 -> 1x4
  _t0_58 = _mm256_blend_pd(_t0_8, _mm256_permute2f128_pd(_mm256_blend_pd(_mm256_setzero_pd(), _t0_7, 1), _mm256_blend_pd(_mm256_setzero_pd(), _t0_7, 1), 8), 12);

  // 4-BLAC: 1x4 Kro 1x4
  _t0_12 = _mm256_mul_pd(_t0_57, _t0_58);

  // 4-BLAC: -( 1x4 )
  _t0_13 = _mm256_sub_pd(_mm256_setzero_pd(), _t0_12);
  _t0_11 = _t0_13;

  // 4x4 -> 4x4 - LowTriang
  _t0_14 = _t0_0;
  _t0_15 = _mm256_blend_pd(_mm256_unpacklo_pd(_t0_5, _t0_2), _mm256_setzero_pd(), 12);
  _t0_16 = _mm256_blend_pd(_t0_9, _mm256_permute2f128_pd(_mm256_blend_pd(_mm256_setzero_pd(), _t0_6, 1), _mm256_blend_pd(_mm256_setzero_pd(), _t0_6, 1), 8), 12);
  _t0_17 = _mm256_blend_pd(_t0_11, _mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_10, _t0_10, 8), _mm256_permute2f128_pd(_t0_10, _t0_10, 8), 8), 8);


  for( int i48 = 0; i48 <= 23; i48+=4 ) {
    _t1_15 = _mm256_broadcast_sd(L + 28*i48 + 112);
    _t1_14 = _mm256_broadcast_sd(L + 28*i48 + 113);
    _t1_13 = _mm256_broadcast_sd(L + 28*i48 + 114);
    _t1_12 = _mm256_broadcast_sd(L + 28*i48 + 115);
    _t1_11 = _mm256_broadcast_sd(L + 28*i48 + 140);
    _t1_10 = _mm256_broadcast_sd(L + 28*i48 + 141);
    _t1_9 = _mm256_broadcast_sd(L + 28*i48 + 142);
    _t1_8 = _mm256_broadcast_sd(L + 28*i48 + 143);
    _t1_7 = _mm256_broadcast_sd(L + 28*i48 + 168);
    _t1_6 = _mm256_broadcast_sd(L + 28*i48 + 169);
    _t1_5 = _mm256_broadcast_sd(L + 28*i48 + 170);
    _t1_4 = _mm256_broadcast_sd(L + 28*i48 + 171);
    _t1_3 = _mm256_broadcast_sd(L + 28*i48 + 196);
    _t1_2 = _mm256_broadcast_sd(L + 28*i48 + 197);
    _t1_1 = _mm256_broadcast_sd(L + 28*i48 + 198);
    _t1_0 = _mm256_broadcast_sd(L + 28*i48 + 199);

    // 4-BLAC: 4x4 * 4x4
    _t1_16 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t1_15, _t0_14), _mm256_mul_pd(_t1_14, _t0_15)), _mm256_add_pd(_mm256_mul_pd(_t1_13, _t0_16), _mm256_mul_pd(_t1_12, _t0_17)));
    _t1_17 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t1_11, _t0_14), _mm256_mul_pd(_t1_10, _t0_15)), _mm256_add_pd(_mm256_mul_pd(_t1_9, _t0_16), _mm256_mul_pd(_t1_8, _t0_17)));
    _t1_18 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t1_7, _t0_14), _mm256_mul_pd(_t1_6, _t0_15)), _mm256_add_pd(_mm256_mul_pd(_t1_5, _t0_16), _mm256_mul_pd(_t1_4, _t0_17)));
    _t1_19 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t1_3, _t0_14), _mm256_mul_pd(_t1_2, _t0_15)), _mm256_add_pd(_mm256_mul_pd(_t1_1, _t0_16), _mm256_mul_pd(_t1_0, _t0_17)));
    _mm256_storeu_pd(L + 28*i48 + 112, _t1_16);
    _mm256_storeu_pd(L + 28*i48 + 140, _t1_17);
    _mm256_storeu_pd(L + 28*i48 + 168, _t1_18);
    _mm256_storeu_pd(L + 28*i48 + 196, _t1_19);
  }


  for( int fi210 = 4; fi210 <= 23; fi210+=4 ) {
    _t2_0 = _mm256_maskload_pd(L + 29*fi210, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
    _t2_1 = _mm256_permute2f128_pd(_mm256_unpacklo_pd(_mm256_maskload_pd(L + 29*fi210 + 28, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0)), _mm256_maskload_pd(L + 29*fi210 + 56, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0))), _mm256_maskload_pd(L + 29*fi210 + 84, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0)), 32);
    _t2_2 = _mm256_maskload_pd(L + 29*fi210 + 29, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
    _t2_3 = _mm256_shuffle_pd(_mm256_maskload_pd(L + 29*fi210 + 57, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0)), _mm256_maskload_pd(L + 29*fi210 + 85, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0)), 0);
    _t2_6 = _mm256_maskload_pd(L + 29*fi210 + 58, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
    _t2_7 = _mm256_maskload_pd(L + 29*fi210 + 86, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
    _t2_10 = _mm256_maskload_pd(L + 29*fi210 + 87, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));

    // Constant 1x1 -> 1x4
    _t2_12 = _mm256_set_pd(0, 0, 0, 1);

    // 1x1 -> 1x4
    _t2_13 = _t2_0;

    // 4-BLAC: 1x4 / 1x4
    _t2_14 = _mm256_div_pd(_t2_12, _t2_13);
    _t2_0 = _t2_14;

    // 3x1 -> 4x1
    _t2_15 = _t2_1;

    // 1x1 -> 1x4
    _t2_16 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t2_0, _t2_0, 32), _mm256_permute2f128_pd(_t2_0, _t2_0, 32), 0);

    // 4-BLAC: 4x1 Kro 1x4
    _t2_17 = _mm256_mul_pd(_t2_15, _t2_16);
    _t2_1 = _t2_17;

    // Constant 1x1 -> 1x4
    _t2_18 = _mm256_set_pd(0, 0, 0, 1);

    // 1x1 -> 1x4
    _t2_19 = _t2_2;

    // 4-BLAC: 1x4 / 1x4
    _t2_20 = _mm256_div_pd(_t2_18, _t2_19);
    _t2_2 = _t2_20;

    // 2x1 -> 4x1
    _t2_21 = _t2_3;

    // 1x1 -> 1x4
    _t2_22 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t2_2, _t2_2, 32), _mm256_permute2f128_pd(_t2_2, _t2_2, 32), 0);

    // 4-BLAC: 4x1 Kro 1x4
    _t2_23 = _mm256_mul_pd(_t2_21, _t2_22);
    _t2_3 = _t2_23;

    // 2x1 -> 4x1
    _t2_24 = _mm256_shuffle_pd(_mm256_blend_pd(_mm256_setzero_pd(), _t2_1, 2), _mm256_permute2f128_pd(_t2_1, _t2_1, 129), 5);

    // 2x1 -> 4x1
    _t2_25 = _t2_3;

    // 1x1 -> 1x4
    _t2_26 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t2_1, _t2_1, 32), _mm256_permute2f128_pd(_t2_1, _t2_1, 32), 0);

    // 4-BLAC: 4x1 Kro 1x4
    _t2_27 = _mm256_mul_pd(_t2_25, _t2_26);

    // 4-BLAC: 4x1 - 4x1
    _t2_28 = _mm256_sub_pd(_t2_24, _t2_27);
    _t2_4 = _t2_28;

    // 1x1 -> 1x4
    _t2_29 = _t2_2;

    // 1x1 -> 1x4
    _t2_30 = _mm256_blend_pd(_mm256_setzero_pd(), _t2_1, 1);

    // 4-BLAC: 1x4 Kro 1x4
    _t2_31 = _mm256_mul_pd(_t2_29, _t2_30);

    // 4-BLAC: -( 1x4 )
    _t2_32 = _mm256_sub_pd(_mm256_setzero_pd(), _t2_31);
    _t2_5 = _t2_32;

    // Constant 1x1 -> 1x4
    _t2_33 = _mm256_set_pd(0, 0, 0, 1);

    // 1x1 -> 1x4
    _t2_34 = _t2_6;

    // 4-BLAC: 1x4 / 1x4
    _t2_35 = _mm256_div_pd(_t2_33, _t2_34);
    _t2_6 = _t2_35;

    // 1x1 -> 1x4
    _t2_36 = _t2_7;

    // 1x1 -> 1x4
    _t2_37 = _t2_6;

    // 4-BLAC: 1x4 Kro 1x4
    _t2_38 = _mm256_mul_pd(_t2_36, _t2_37);
    _t2_7 = _t2_38;

    // 1x2 -> 1x4
    _t2_39 = _mm256_unpackhi_pd(_mm256_blend_pd(_t2_4, _mm256_setzero_pd(), 12), _mm256_blend_pd(_t2_3, _mm256_setzero_pd(), 12));

    // 1x1 -> 1x4
    _t2_40 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t2_7, _t2_7, 32), _mm256_permute2f128_pd(_t2_7, _t2_7, 32), 0);

    // 1x2 -> 1x4
    _t2_41 = _mm256_blend_pd(_mm256_unpacklo_pd(_t2_4, _t2_3), _mm256_setzero_pd(), 12);

    // 4-BLAC: 1x4 Kro 1x4
    _t2_42 = _mm256_mul_pd(_t2_40, _t2_41);

    // 4-BLAC: 1x4 - 1x4
    _t2_43 = _mm256_sub_pd(_t2_39, _t2_42);
    _t2_8 = _t2_43;

    // 1x1 -> 1x4
    _t2_44 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t2_6, _t2_6, 32), _mm256_permute2f128_pd(_t2_6, _t2_6, 32), 0);

    // 1x2 -> 1x4
    _t2_45 = _mm256_blend_pd(_mm256_unpacklo_pd(_t2_4, _t2_3), _mm256_setzero_pd(), 12);

    // 4-BLAC: 1x4 Kro 1x4
    _t2_46 = _mm256_mul_pd(_t2_44, _t2_45);

    // 4-BLAC: -( 1x4 )
    _t2_47 = _mm256_sub_pd(_mm256_setzero_pd(), _t2_46);
    _t2_9 = _t2_47;

    // Constant 1x1 -> 1x4
    _t2_48 = _mm256_set_pd(0, 0, 0, 1);

    // 1x1 -> 1x4
    _t2_49 = _t2_10;

    // 4-BLAC: 1x4 / 1x4
    _t2_50 = _mm256_div_pd(_t2_48, _t2_49);
    _t2_10 = _t2_50;

    // 1x1 -> 1x4
    _t2_51 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t2_10, _t2_10, 32), _mm256_permute2f128_pd(_t2_10, _t2_10, 32), 0);

    // 1x3 -> 1x4
    _t2_52 = _mm256_blend_pd(_t2_8, _mm256_permute2f128_pd(_mm256_blend_pd(_mm256_setzero_pd(), _t2_7, 1), _mm256_blend_pd(_mm256_setzero_pd(), _t2_7, 1), 8), 12);

    // 4-BLAC: 1x4 Kro 1x4
    _t2_53 = _mm256_mul_pd(_t2_51, _t2_52);

    // 4-BLAC: -( 1x4 )
    _t2_54 = _mm256_sub_pd(_mm256_setzero_pd(), _t2_53);
    _t2_11 = _t2_54;

    // 4x4 -> 4x4 - LowTriang
    _t2_55 = _t2_0;
    _t2_56 = _mm256_blend_pd(_mm256_unpacklo_pd(_t2_5, _t2_2), _mm256_setzero_pd(), 12);
    _t2_57 = _mm256_blend_pd(_t2_9, _mm256_permute2f128_pd(_mm256_blend_pd(_mm256_setzero_pd(), _t2_6, 1), _mm256_blend_pd(_mm256_setzero_pd(), _t2_6, 1), 8), 12);
    _t2_58 = _mm256_blend_pd(_t2_11, _mm256_shuffle_pd(_mm256_permute2f128_pd(_t2_10, _t2_10, 8), _mm256_permute2f128_pd(_t2_10, _t2_10, 8), 8), 8);
    _mm256_maskstore_pd(L + 29*fi210 + 28, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_1);
    _mm256_maskstore_pd(L + 29*fi210 + 56, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t2_1, _t2_1, 1));
    _mm256_maskstore_pd(L + 29*fi210 + 84, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_permute2f128_pd(_t2_1, _t2_1, 129));
    _mm256_maskstore_pd(L + 29*fi210 + 57, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_3);
    _mm256_maskstore_pd(L + 29*fi210 + 85, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t2_3, _t2_3, 1));
    _mm256_maskstore_pd(L + 29*fi210 + 56, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_4);
    _mm256_maskstore_pd(L + 29*fi210 + 84, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t2_4, _t2_4, 1));
    _mm256_maskstore_pd(L + 29*fi210 + 86, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_7);
    _mm256_maskstore_pd(L + 29*fi210 + 84, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, 0, 0), _t2_8);

    for( int i48 = 0; i48 <= -fi210 + 23; i48+=4 ) {
      _t3_15 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 112);
      _t3_14 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 113);
      _t3_13 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 114);
      _t3_12 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 115);
      _t3_11 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 140);
      _t3_10 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 141);
      _t3_9 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 142);
      _t3_8 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 143);
      _t3_7 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 168);
      _t3_6 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 169);
      _t3_5 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 170);
      _t3_4 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 171);
      _t3_3 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 196);
      _t3_2 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 197);
      _t3_1 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 198);
      _t3_0 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 199);

      // 4x4 -> 4x4 - LowTriang
      _t3_20 = _t2_0;
      _t3_21 = _mm256_blend_pd(_mm256_unpacklo_pd(_t2_5, _t2_2), _mm256_setzero_pd(), 12);
      _t3_22 = _mm256_blend_pd(_t2_9, _mm256_permute2f128_pd(_mm256_blend_pd(_mm256_setzero_pd(), _t2_6, 1), _mm256_blend_pd(_mm256_setzero_pd(), _t2_6, 1), 8), 12);
      _t3_23 = _mm256_blend_pd(_t2_11, _mm256_shuffle_pd(_mm256_permute2f128_pd(_t2_10, _t2_10, 8), _mm256_permute2f128_pd(_t2_10, _t2_10, 8), 8), 8);

      // 4-BLAC: 4x4 * 4x4
      _t3_16 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t3_15, _t3_20), _mm256_mul_pd(_t3_14, _t3_21)), _mm256_add_pd(_mm256_mul_pd(_t3_13, _t3_22), _mm256_mul_pd(_t3_12, _t3_23)));
      _t3_17 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t3_11, _t3_20), _mm256_mul_pd(_t3_10, _t3_21)), _mm256_add_pd(_mm256_mul_pd(_t3_9, _t3_22), _mm256_mul_pd(_t3_8, _t3_23)));
      _t3_18 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t3_7, _t3_20), _mm256_mul_pd(_t3_6, _t3_21)), _mm256_add_pd(_mm256_mul_pd(_t3_5, _t3_22), _mm256_mul_pd(_t3_4, _t3_23)));
      _t3_19 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t3_3, _t3_20), _mm256_mul_pd(_t3_2, _t3_21)), _mm256_add_pd(_mm256_mul_pd(_t3_1, _t3_22), _mm256_mul_pd(_t3_0, _t3_23)));
      _mm256_storeu_pd(L + 29*fi210 + 28*i48 + 112, _t3_16);
      _mm256_storeu_pd(L + 29*fi210 + 28*i48 + 140, _t3_17);
      _mm256_storeu_pd(L + 29*fi210 + 28*i48 + 168, _t3_18);
      _mm256_storeu_pd(L + 29*fi210 + 28*i48 + 196, _t3_19);
    }

    for( int i48 = 0; i48 <= -fi210 + 23; i48+=4 ) {

      for( int i99 = 0; i99 <= fi210 - 1; i99+=4 ) {
        _t4_24 = _mm256_loadu_pd(L + 28*fi210 + 28*i48 + i99 + 112);
        _t4_25 = _mm256_loadu_pd(L + 28*fi210 + 28*i48 + i99 + 140);
        _t4_26 = _mm256_loadu_pd(L + 28*fi210 + 28*i48 + i99 + 168);
        _t4_27 = _mm256_loadu_pd(L + 28*fi210 + 28*i48 + i99 + 196);
        _t4_19 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 112);
        _t4_18 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 113);
        _t4_17 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 114);
        _t4_16 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 115);
        _t4_15 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 140);
        _t4_14 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 141);
        _t4_13 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 142);
        _t4_12 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 143);
        _t4_11 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 168);
        _t4_10 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 169);
        _t4_9 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 170);
        _t4_8 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 171);
        _t4_7 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 196);
        _t4_6 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 197);
        _t4_5 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 198);
        _t4_4 = _mm256_broadcast_sd(L + 29*fi210 + 28*i48 + 199);
        _t4_3 = _mm256_loadu_pd(L + 28*fi210 + i99);
        _t4_2 = _mm256_loadu_pd(L + 28*fi210 + i99 + 28);
        _t4_1 = _mm256_loadu_pd(L + 28*fi210 + i99 + 56);
        _t4_0 = _mm256_loadu_pd(L + 28*fi210 + i99 + 84);

        // 4-BLAC: 4x4 * 4x4
        _t4_20 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t4_19, _t4_3), _mm256_mul_pd(_t4_18, _t4_2)), _mm256_add_pd(_mm256_mul_pd(_t4_17, _t4_1), _mm256_mul_pd(_t4_16, _t4_0)));
        _t4_21 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t4_15, _t4_3), _mm256_mul_pd(_t4_14, _t4_2)), _mm256_add_pd(_mm256_mul_pd(_t4_13, _t4_1), _mm256_mul_pd(_t4_12, _t4_0)));
        _t4_22 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t4_11, _t4_3), _mm256_mul_pd(_t4_10, _t4_2)), _mm256_add_pd(_mm256_mul_pd(_t4_9, _t4_1), _mm256_mul_pd(_t4_8, _t4_0)));
        _t4_23 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t4_7, _t4_3), _mm256_mul_pd(_t4_6, _t4_2)), _mm256_add_pd(_mm256_mul_pd(_t4_5, _t4_1), _mm256_mul_pd(_t4_4, _t4_0)));

        // 4-BLAC: 4x4 - 4x4
        _t4_24 = _mm256_sub_pd(_t4_24, _t4_20);
        _t4_25 = _mm256_sub_pd(_t4_25, _t4_21);
        _t4_26 = _mm256_sub_pd(_t4_26, _t4_22);
        _t4_27 = _mm256_sub_pd(_t4_27, _t4_23);
        _mm256_storeu_pd(L + 28*fi210 + 28*i48 + i99 + 112, _t4_24);
        _mm256_storeu_pd(L + 28*fi210 + 28*i48 + i99 + 140, _t4_25);
        _mm256_storeu_pd(L + 28*fi210 + 28*i48 + i99 + 168, _t4_26);
        _mm256_storeu_pd(L + 28*fi210 + 28*i48 + i99 + 196, _t4_27);
      }
    }

    // 4x4 -> 4x4 - LowTriang
    _t5_0 = _t2_0;
    _t5_1 = _mm256_blend_pd(_mm256_unpacklo_pd(_t2_5, _t2_2), _mm256_setzero_pd(), 12);
    _t5_2 = _mm256_blend_pd(_t2_9, _mm256_permute2f128_pd(_mm256_blend_pd(_mm256_setzero_pd(), _t2_6, 1), _mm256_blend_pd(_mm256_setzero_pd(), _t2_6, 1), 8), 12);
    _t5_3 = _mm256_blend_pd(_t2_11, _mm256_shuffle_pd(_mm256_permute2f128_pd(_t2_10, _t2_10, 8), _mm256_permute2f128_pd(_t2_10, _t2_10, 8), 8), 8);

    for( int i48 = 0; i48 <= fi210 - 1; i48+=4 ) {
      _t6_4 = _mm256_loadu_pd(L + 28*fi210 + i48);
      _t6_5 = _mm256_loadu_pd(L + 28*fi210 + i48 + 28);
      _t6_6 = _mm256_loadu_pd(L + 28*fi210 + i48 + 56);
      _t6_7 = _mm256_loadu_pd(L + 28*fi210 + i48 + 84);

      // 4x4 -> 4x4 - LowTriang
      _t6_8 = _t2_0;
      _t6_9 = _mm256_blend_pd(_mm256_unpacklo_pd(_t2_5, _t2_2), _mm256_setzero_pd(), 12);
      _t6_10 = _mm256_blend_pd(_t2_9, _mm256_permute2f128_pd(_mm256_blend_pd(_mm256_setzero_pd(), _t2_6, 1), _mm256_blend_pd(_mm256_setzero_pd(), _t2_6, 1), 8), 12);
      _t6_11 = _mm256_blend_pd(_t2_11, _mm256_shuffle_pd(_mm256_permute2f128_pd(_t2_10, _t2_10, 8), _mm256_permute2f128_pd(_t2_10, _t2_10, 8), 8), 8);

      // 4-BLAC: 4x4 * 4x4
      _t6_0 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t6_8, _t6_8, 32), _mm256_permute2f128_pd(_t6_8, _t6_8, 32), 0), _t6_4), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t6_8, _t6_8, 32), _mm256_permute2f128_pd(_t6_8, _t6_8, 32), 15), _t6_5)), _mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t6_8, _t6_8, 49), _mm256_permute2f128_pd(_t6_8, _t6_8, 49), 0), _t6_6), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t6_8, _t6_8, 49), _mm256_permute2f128_pd(_t6_8, _t6_8, 49), 15), _t6_7)));
      _t6_1 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t6_9, _t6_9, 32), _mm256_permute2f128_pd(_t6_9, _t6_9, 32), 0), _t6_4), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t6_9, _t6_9, 32), _mm256_permute2f128_pd(_t6_9, _t6_9, 32), 15), _t6_5)), _mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t6_9, _t6_9, 49), _mm256_permute2f128_pd(_t6_9, _t6_9, 49), 0), _t6_6), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t6_9, _t6_9, 49), _mm256_permute2f128_pd(_t6_9, _t6_9, 49), 15), _t6_7)));
      _t6_2 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t6_10, _t6_10, 32), _mm256_permute2f128_pd(_t6_10, _t6_10, 32), 0), _t6_4), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t6_10, _t6_10, 32), _mm256_permute2f128_pd(_t6_10, _t6_10, 32), 15), _t6_5)), _mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t6_10, _t6_10, 49), _mm256_permute2f128_pd(_t6_10, _t6_10, 49), 0), _t6_6), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t6_10, _t6_10, 49), _mm256_permute2f128_pd(_t6_10, _t6_10, 49), 15), _t6_7)));
      _t6_3 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t6_11, _t6_11, 32), _mm256_permute2f128_pd(_t6_11, _t6_11, 32), 0), _t6_4), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t6_11, _t6_11, 32), _mm256_permute2f128_pd(_t6_11, _t6_11, 32), 15), _t6_5)), _mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t6_11, _t6_11, 49), _mm256_permute2f128_pd(_t6_11, _t6_11, 49), 0), _t6_6), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t6_11, _t6_11, 49), _mm256_permute2f128_pd(_t6_11, _t6_11, 49), 15), _t6_7)));

      // 4-BLAC: -( 4x4 )
      _t6_4 = _mm256_sub_pd(_mm256_setzero_pd(), _t6_0);
      _t6_5 = _mm256_sub_pd(_mm256_setzero_pd(), _t6_1);
      _t6_6 = _mm256_sub_pd(_mm256_setzero_pd(), _t6_2);
      _t6_7 = _mm256_sub_pd(_mm256_setzero_pd(), _t6_3);
      _mm256_storeu_pd(L + 28*fi210 + i48, _t6_4);
      _mm256_storeu_pd(L + 28*fi210 + i48 + 28, _t6_5);
      _mm256_storeu_pd(L + 28*fi210 + i48 + 56, _t6_6);
      _mm256_storeu_pd(L + 28*fi210 + i48 + 84, _t6_7);
    }
    _mm256_maskstore_pd(L + 29*fi210, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_0);
    _mm256_maskstore_pd(L + 29*fi210 + 29, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_2);
    _mm256_maskstore_pd(L + 29*fi210 + 28, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_5);
    _mm256_maskstore_pd(L + 29*fi210 + 58, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_6);
    _mm256_maskstore_pd(L + 29*fi210 + 56, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, 0, 0), _t2_9);
    _mm256_maskstore_pd(L + 29*fi210 + 87, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_10);
    _mm256_maskstore_pd(L + 29*fi210 + 84, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, (__int64)1 << 63, 0), _t2_11);
  }

  _t7_0 = _mm256_maskload_pd(L + 696, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
  _t7_1 = _mm256_permute2f128_pd(_mm256_unpacklo_pd(_mm256_maskload_pd(L + 724, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0)), _mm256_maskload_pd(L + 752, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0))), _mm256_maskload_pd(L + 780, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0)), 32);
  _t7_2 = _mm256_maskload_pd(L + 725, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
  _t7_3 = _mm256_shuffle_pd(_mm256_maskload_pd(L + 753, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0)), _mm256_maskload_pd(L + 781, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0)), 0);
  _t7_6 = _mm256_maskload_pd(L + 754, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
  _t7_7 = _mm256_maskload_pd(L + 782, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
  _t7_10 = _mm256_maskload_pd(L + 783, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));

  // Constant 1x1 -> 1x4
  _t7_12 = _mm256_set_pd(0, 0, 0, 1);

  // 1x1 -> 1x4
  _t7_13 = _t7_0;

  // 4-BLAC: 1x4 / 1x4
  _t7_14 = _mm256_div_pd(_t7_12, _t7_13);
  _t7_0 = _t7_14;

  // 3x1 -> 4x1
  _t7_15 = _t7_1;

  // 1x1 -> 1x4
  _t7_16 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_0, _t7_0, 32), _mm256_permute2f128_pd(_t7_0, _t7_0, 32), 0);

  // 4-BLAC: 4x1 Kro 1x4
  _t7_17 = _mm256_mul_pd(_t7_15, _t7_16);
  _t7_1 = _t7_17;

  // Constant 1x1 -> 1x4
  _t7_18 = _mm256_set_pd(0, 0, 0, 1);

  // 1x1 -> 1x4
  _t7_19 = _t7_2;

  // 4-BLAC: 1x4 / 1x4
  _t7_20 = _mm256_div_pd(_t7_18, _t7_19);
  _t7_2 = _t7_20;

  // 2x1 -> 4x1
  _t7_21 = _t7_3;

  // 1x1 -> 1x4
  _t7_22 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_2, _t7_2, 32), _mm256_permute2f128_pd(_t7_2, _t7_2, 32), 0);

  // 4-BLAC: 4x1 Kro 1x4
  _t7_23 = _mm256_mul_pd(_t7_21, _t7_22);
  _t7_3 = _t7_23;

  // 2x1 -> 4x1
  _t7_24 = _mm256_shuffle_pd(_mm256_blend_pd(_mm256_setzero_pd(), _t7_1, 2), _mm256_permute2f128_pd(_t7_1, _t7_1, 129), 5);

  // 2x1 -> 4x1
  _t7_25 = _t7_3;

  // 1x1 -> 1x4
  _t7_26 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_1, _t7_1, 32), _mm256_permute2f128_pd(_t7_1, _t7_1, 32), 0);

  // 4-BLAC: 4x1 Kro 1x4
  _t7_27 = _mm256_mul_pd(_t7_25, _t7_26);

  // 4-BLAC: 4x1 - 4x1
  _t7_28 = _mm256_sub_pd(_t7_24, _t7_27);
  _t7_4 = _t7_28;

  // 1x1 -> 1x4
  _t7_29 = _t7_2;

  // 1x1 -> 1x4
  _t7_30 = _mm256_blend_pd(_mm256_setzero_pd(), _t7_1, 1);

  // 4-BLAC: 1x4 Kro 1x4
  _t7_31 = _mm256_mul_pd(_t7_29, _t7_30);

  // 4-BLAC: -( 1x4 )
  _t7_32 = _mm256_sub_pd(_mm256_setzero_pd(), _t7_31);
  _t7_5 = _t7_32;

  // Constant 1x1 -> 1x4
  _t7_33 = _mm256_set_pd(0, 0, 0, 1);

  // 1x1 -> 1x4
  _t7_34 = _t7_6;

  // 4-BLAC: 1x4 / 1x4
  _t7_35 = _mm256_div_pd(_t7_33, _t7_34);
  _t7_6 = _t7_35;

  // 1x1 -> 1x4
  _t7_36 = _t7_7;

  // 1x1 -> 1x4
  _t7_37 = _t7_6;

  // 4-BLAC: 1x4 Kro 1x4
  _t7_38 = _mm256_mul_pd(_t7_36, _t7_37);
  _t7_7 = _t7_38;

  // 1x2 -> 1x4
  _t7_39 = _mm256_unpackhi_pd(_mm256_blend_pd(_t7_4, _mm256_setzero_pd(), 12), _mm256_blend_pd(_t7_3, _mm256_setzero_pd(), 12));

  // 1x1 -> 1x4
  _t7_40 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_7, _t7_7, 32), _mm256_permute2f128_pd(_t7_7, _t7_7, 32), 0);

  // 1x2 -> 1x4
  _t7_41 = _mm256_blend_pd(_mm256_unpacklo_pd(_t7_4, _t7_3), _mm256_setzero_pd(), 12);

  // 4-BLAC: 1x4 Kro 1x4
  _t7_42 = _mm256_mul_pd(_t7_40, _t7_41);

  // 4-BLAC: 1x4 - 1x4
  _t7_43 = _mm256_sub_pd(_t7_39, _t7_42);
  _t7_8 = _t7_43;

  // 1x1 -> 1x4
  _t7_44 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_6, _t7_6, 32), _mm256_permute2f128_pd(_t7_6, _t7_6, 32), 0);

  // 1x2 -> 1x4
  _t7_45 = _mm256_blend_pd(_mm256_unpacklo_pd(_t7_4, _t7_3), _mm256_setzero_pd(), 12);

  // 4-BLAC: 1x4 Kro 1x4
  _t7_46 = _mm256_mul_pd(_t7_44, _t7_45);

  // 4-BLAC: -( 1x4 )
  _t7_47 = _mm256_sub_pd(_mm256_setzero_pd(), _t7_46);
  _t7_9 = _t7_47;

  // Constant 1x1 -> 1x4
  _t7_48 = _mm256_set_pd(0, 0, 0, 1);

  // 1x1 -> 1x4
  _t7_49 = _t7_10;

  // 4-BLAC: 1x4 / 1x4
  _t7_50 = _mm256_div_pd(_t7_48, _t7_49);
  _t7_10 = _t7_50;

  // 1x1 -> 1x4
  _t7_51 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_10, _t7_10, 32), _mm256_permute2f128_pd(_t7_10, _t7_10, 32), 0);

  // 1x3 -> 1x4
  _t7_52 = _mm256_blend_pd(_t7_8, _mm256_permute2f128_pd(_mm256_blend_pd(_mm256_setzero_pd(), _t7_7, 1), _mm256_blend_pd(_mm256_setzero_pd(), _t7_7, 1), 8), 12);

  // 4-BLAC: 1x4 Kro 1x4
  _t7_53 = _mm256_mul_pd(_t7_51, _t7_52);

  // 4-BLAC: -( 1x4 )
  _t7_54 = _mm256_sub_pd(_mm256_setzero_pd(), _t7_53);
  _t7_11 = _t7_54;

  // 4x4 -> 4x4 - LowTriang
  _t7_55 = _t7_0;
  _t7_56 = _mm256_blend_pd(_mm256_unpacklo_pd(_t7_5, _t7_2), _mm256_setzero_pd(), 12);
  _t7_57 = _mm256_blend_pd(_t7_9, _mm256_permute2f128_pd(_mm256_blend_pd(_mm256_setzero_pd(), _t7_6, 1), _mm256_blend_pd(_mm256_setzero_pd(), _t7_6, 1), 8), 12);
  _t7_58 = _mm256_blend_pd(_t7_11, _mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_10, _t7_10, 8), _mm256_permute2f128_pd(_t7_10, _t7_10, 8), 8), 8);


  for( int i48 = 0; i48 <= 23; i48+=4 ) {
    _t8_4 = _mm256_loadu_pd(L + i48 + 672);
    _t8_5 = _mm256_loadu_pd(L + i48 + 700);
    _t8_6 = _mm256_loadu_pd(L + i48 + 728);
    _t8_7 = _mm256_loadu_pd(L + i48 + 756);

    // 4-BLAC: 4x4 * 4x4
    _t8_0 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_55, _t7_55, 32), _mm256_permute2f128_pd(_t7_55, _t7_55, 32), 0), _t8_4), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_55, _t7_55, 32), _mm256_permute2f128_pd(_t7_55, _t7_55, 32), 15), _t8_5)), _mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_55, _t7_55, 49), _mm256_permute2f128_pd(_t7_55, _t7_55, 49), 0), _t8_6), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_55, _t7_55, 49), _mm256_permute2f128_pd(_t7_55, _t7_55, 49), 15), _t8_7)));
    _t8_1 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_56, _t7_56, 32), _mm256_permute2f128_pd(_t7_56, _t7_56, 32), 0), _t8_4), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_56, _t7_56, 32), _mm256_permute2f128_pd(_t7_56, _t7_56, 32), 15), _t8_5)), _mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_56, _t7_56, 49), _mm256_permute2f128_pd(_t7_56, _t7_56, 49), 0), _t8_6), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_56, _t7_56, 49), _mm256_permute2f128_pd(_t7_56, _t7_56, 49), 15), _t8_7)));
    _t8_2 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_57, _t7_57, 32), _mm256_permute2f128_pd(_t7_57, _t7_57, 32), 0), _t8_4), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_57, _t7_57, 32), _mm256_permute2f128_pd(_t7_57, _t7_57, 32), 15), _t8_5)), _mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_57, _t7_57, 49), _mm256_permute2f128_pd(_t7_57, _t7_57, 49), 0), _t8_6), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_57, _t7_57, 49), _mm256_permute2f128_pd(_t7_57, _t7_57, 49), 15), _t8_7)));
    _t8_3 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_58, _t7_58, 32), _mm256_permute2f128_pd(_t7_58, _t7_58, 32), 0), _t8_4), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_58, _t7_58, 32), _mm256_permute2f128_pd(_t7_58, _t7_58, 32), 15), _t8_5)), _mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_58, _t7_58, 49), _mm256_permute2f128_pd(_t7_58, _t7_58, 49), 0), _t8_6), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_58, _t7_58, 49), _mm256_permute2f128_pd(_t7_58, _t7_58, 49), 15), _t8_7)));

    // 4-BLAC: -( 4x4 )
    _t8_4 = _mm256_sub_pd(_mm256_setzero_pd(), _t8_0);
    _t8_5 = _mm256_sub_pd(_mm256_setzero_pd(), _t8_1);
    _t8_6 = _mm256_sub_pd(_mm256_setzero_pd(), _t8_2);
    _t8_7 = _mm256_sub_pd(_mm256_setzero_pd(), _t8_3);
    _mm256_storeu_pd(L + i48 + 672, _t8_4);
    _mm256_storeu_pd(L + i48 + 700, _t8_5);
    _mm256_storeu_pd(L + i48 + 728, _t8_6);
    _mm256_storeu_pd(L + i48 + 756, _t8_7);
  }

  _mm256_maskstore_pd(L, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_0);
  _mm256_maskstore_pd(L + 28, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_1);
  _mm256_maskstore_pd(L + 56, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t0_1, _t0_1, 1));
  _mm256_maskstore_pd(L + 84, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_permute2f128_pd(_t0_1, _t0_1, 129));
  _mm256_maskstore_pd(L + 29, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_2);
  _mm256_maskstore_pd(L + 57, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_3);
  _mm256_maskstore_pd(L + 85, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t0_3, _t0_3, 1));
  _mm256_maskstore_pd(L + 56, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_4);
  _mm256_maskstore_pd(L + 84, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t0_4, _t0_4, 1));
  _mm256_maskstore_pd(L + 28, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_5);
  _mm256_maskstore_pd(L + 58, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_6);
  _mm256_maskstore_pd(L + 86, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_7);
  _mm256_maskstore_pd(L + 84, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, 0, 0), _t0_8);
  _mm256_maskstore_pd(L + 56, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, 0, 0), _t0_9);
  _mm256_maskstore_pd(L + 87, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_10);
  _mm256_maskstore_pd(L + 84, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, (__int64)1 << 63, 0), _t0_11);
  _mm256_maskstore_pd(L + 696, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_0);
  _mm256_maskstore_pd(L + 724, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_1);
  _mm256_maskstore_pd(L + 752, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t7_1, _t7_1, 1));
  _mm256_maskstore_pd(L + 780, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_permute2f128_pd(_t7_1, _t7_1, 129));
  _mm256_maskstore_pd(L + 725, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_2);
  _mm256_maskstore_pd(L + 753, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_3);
  _mm256_maskstore_pd(L + 781, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t7_3, _t7_3, 1));
  _mm256_maskstore_pd(L + 752, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_4);
  _mm256_maskstore_pd(L + 780, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t7_4, _t7_4, 1));
  _mm256_maskstore_pd(L + 724, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_5);
  _mm256_maskstore_pd(L + 754, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_6);
  _mm256_maskstore_pd(L + 782, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_7);
  _mm256_maskstore_pd(L + 780, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, 0, 0), _t7_8);
  _mm256_maskstore_pd(L + 752, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, 0, 0), _t7_9);
  _mm256_maskstore_pd(L + 783, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_10);
  _mm256_maskstore_pd(L + 780, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, (__int64)1 << 63, 0), _t7_11);

}
