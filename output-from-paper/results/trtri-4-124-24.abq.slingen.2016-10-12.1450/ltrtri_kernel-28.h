/*
 * ltrtri_kernel.h
 *
Decl { {u'L': LowerTriangular[L, (28, 28), GenMatAccess]} }
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ann: {'part_schemes': {'tril_inv_ow_opt': {'m': 'm4.ll'}}, 'cl1ck_v': 0, 'variant_tag': 'tril_inv_ow_opt_m4'}

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Entry 0:
Eq: Tile( (1, 1), G(h(1, 28, 27), L[28,28],h(1, 28, 27)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, 27), L[28,28],h(1, 28, 27)) ) )
Eq.ann: {}
Entry 1:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 27), L[28,28],h(3, 28, 24)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 27), L[28,28],h(1, 28, 27)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 27), L[28,28],h(3, 28, 24)) ) ) ) )
Eq.ann: {}
Entry 2:
Eq: Tile( (1, 1), G(h(1, 28, 26), L[28,28],h(1, 28, 26)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, 26), L[28,28],h(1, 28, 26)) ) )
Eq.ann: {}
Entry 3:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 27), L[28,28],h(1, 28, 26)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 27), L[28,28],h(1, 28, 26)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 26), L[28,28],h(1, 28, 26)) ) ) )
Eq.ann: {}
Entry 4:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 27), L[28,28],h(2, 28, 24)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 27), L[28,28],h(2, 28, 24)) ) ) - ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 27), L[28,28],h(1, 28, 26)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 26), L[28,28],h(2, 28, 24)) ) ) ) )
Eq.ann: {}
Entry 5:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 26), L[28,28],h(2, 28, 24)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 26), L[28,28],h(1, 28, 26)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 26), L[28,28],h(2, 28, 24)) ) ) ) )
Eq.ann: {}
Entry 6:
Eq: Tile( (1, 1), G(h(1, 28, 25), L[28,28],h(1, 28, 25)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, 25), L[28,28],h(1, 28, 25)) ) )
Eq.ann: {}
Entry 7:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 26), L[28,28],h(1, 28, 25)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 26), L[28,28],h(1, 28, 25)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 25), L[28,28],h(1, 28, 25)) ) ) )
Eq.ann: {}
Entry 8:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 26), L[28,28],h(1, 28, 24)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 26), L[28,28],h(1, 28, 24)) ) ) - ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 26), L[28,28],h(1, 28, 25)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 25), L[28,28],h(1, 28, 24)) ) ) ) )
Eq.ann: {}
Entry 9:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 25), L[28,28],h(1, 28, 24)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 25), L[28,28],h(1, 28, 25)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 25), L[28,28],h(1, 28, 24)) ) ) ) )
Eq.ann: {}
Entry 10:
Eq: Tile( (1, 1), G(h(1, 28, 24), L[28,28],h(1, 28, 24)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, 24), L[28,28],h(1, 28, 24)) ) )
Eq.ann: {}
Entry 11:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(3, 28, 25), L[28,28],h(1, 28, 24)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(3, 28, 25), L[28,28],h(1, 28, 24)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 24), L[28,28],h(1, 28, 24)) ) ) )
Eq.ann: {}
Entry 12:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(4, 28, 24), L[28,28],h(24, 28, 0)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(4, 28, 24), L[28,28],h(4, 28, 24)) ) ) * Tile( (1, 1), Tile( (4, 4), G(h(4, 28, 24), L[28,28],h(24, 28, 0)) ) ) ) )
Eq.ann: {}
Entry 13:
For_{fi2;4;23;4} ( Entry 0:
Eq: Tile( (1, 1), G(h(1, 28, -fi2 + 27), L[28,28],h(1, 28, -fi2 + 27)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, -fi2 + 27), L[28,28],h(1, 28, -fi2 + 27)) ) )
Eq.ann: {}
Entry 1:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 27), L[28,28],h(3, 28, -fi2 + 24)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 27), L[28,28],h(1, 28, -fi2 + 27)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 27), L[28,28],h(3, 28, -fi2 + 24)) ) ) ) )
Eq.ann: {}
Entry 2:
Eq: Tile( (1, 1), G(h(1, 28, -fi2 + 26), L[28,28],h(1, 28, -fi2 + 26)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, -fi2 + 26), L[28,28],h(1, 28, -fi2 + 26)) ) )
Eq.ann: {}
Entry 3:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 27), L[28,28],h(1, 28, -fi2 + 26)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 27), L[28,28],h(1, 28, -fi2 + 26)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 26), L[28,28],h(1, 28, -fi2 + 26)) ) ) )
Eq.ann: {}
Entry 4:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 27), L[28,28],h(2, 28, -fi2 + 24)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 27), L[28,28],h(2, 28, -fi2 + 24)) ) ) - ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 27), L[28,28],h(1, 28, -fi2 + 26)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 26), L[28,28],h(2, 28, -fi2 + 24)) ) ) ) )
Eq.ann: {}
Entry 5:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 26), L[28,28],h(2, 28, -fi2 + 24)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 26), L[28,28],h(1, 28, -fi2 + 26)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 26), L[28,28],h(2, 28, -fi2 + 24)) ) ) ) )
Eq.ann: {}
Entry 6:
Eq: Tile( (1, 1), G(h(1, 28, -fi2 + 25), L[28,28],h(1, 28, -fi2 + 25)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, -fi2 + 25), L[28,28],h(1, 28, -fi2 + 25)) ) )
Eq.ann: {}
Entry 7:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(2, 28, -fi2 + 26), L[28,28],h(1, 28, -fi2 + 25)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, -fi2 + 26), L[28,28],h(1, 28, -fi2 + 25)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 25), L[28,28],h(1, 28, -fi2 + 25)) ) ) )
Eq.ann: {}
Entry 8:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(2, 28, -fi2 + 26), L[28,28],h(1, 28, -fi2 + 24)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, -fi2 + 26), L[28,28],h(1, 28, -fi2 + 24)) ) ) - ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, -fi2 + 26), L[28,28],h(1, 28, -fi2 + 25)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 25), L[28,28],h(1, 28, -fi2 + 24)) ) ) ) )
Eq.ann: {}
Entry 9:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 25), L[28,28],h(1, 28, -fi2 + 24)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 25), L[28,28],h(1, 28, -fi2 + 25)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 25), L[28,28],h(1, 28, -fi2 + 24)) ) ) ) )
Eq.ann: {}
Entry 10:
Eq: Tile( (1, 1), G(h(1, 28, -fi2 + 24), L[28,28],h(1, 28, -fi2 + 24)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, -fi2 + 24), L[28,28],h(1, 28, -fi2 + 24)) ) )
Eq.ann: {}
Entry 11:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(3, 28, -fi2 + 25), L[28,28],h(1, 28, -fi2 + 24)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(3, 28, -fi2 + 25), L[28,28],h(1, 28, -fi2 + 24)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, -fi2 + 24), L[28,28],h(1, 28, -fi2 + 24)) ) ) )
Eq.ann: {}
Entry 12:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(fi2, 28, -fi2 + 28), L[28,28],h(4, 28, -fi2 + 24)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(fi2, 28, -fi2 + 28), L[28,28],h(4, 28, -fi2 + 24)) ) ) * Tile( (1, 1), Tile( (4, 4), G(h(4, 28, -fi2 + 24), L[28,28],h(4, 28, -fi2 + 24)) ) ) )
Eq.ann: {}
Entry 13:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(fi2, 28, -fi2 + 28), L[28,28],h(-fi2 + 24, 28, 0)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(fi2, 28, -fi2 + 28), L[28,28],h(-fi2 + 24, 28, 0)) ) ) - ( Tile( (1, 1), Tile( (4, 4), G(h(fi2, 28, -fi2 + 28), L[28,28],h(4, 28, -fi2 + 24)) ) ) * Tile( (1, 1), Tile( (4, 4), G(h(4, 28, -fi2 + 24), L[28,28],h(-fi2 + 24, 28, 0)) ) ) ) )
Eq.ann: {}
Entry 14:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(4, 28, -fi2 + 24), L[28,28],h(-fi2 + 24, 28, 0)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(4, 28, -fi2 + 24), L[28,28],h(4, 28, -fi2 + 24)) ) ) * Tile( (1, 1), Tile( (4, 4), G(h(4, 28, -fi2 + 24), L[28,28],h(-fi2 + 24, 28, 0)) ) ) ) )
Eq.ann: {}
 )Entry 14:
Eq: Tile( (1, 1), G(h(1, 28, 3), L[28,28],h(1, 28, 3)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, 3), L[28,28],h(1, 28, 3)) ) )
Eq.ann: {}
Entry 15:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 3), L[28,28],h(3, 28, 0)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 3), L[28,28],h(1, 28, 3)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 3), L[28,28],h(3, 28, 0)) ) ) ) )
Eq.ann: {}
Entry 16:
Eq: Tile( (1, 1), G(h(1, 28, 2), L[28,28],h(1, 28, 2)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, 2), L[28,28],h(1, 28, 2)) ) )
Eq.ann: {}
Entry 17:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 3), L[28,28],h(1, 28, 2)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 3), L[28,28],h(1, 28, 2)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 2), L[28,28],h(1, 28, 2)) ) ) )
Eq.ann: {}
Entry 18:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 3), L[28,28],h(2, 28, 0)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 3), L[28,28],h(2, 28, 0)) ) ) - ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 3), L[28,28],h(1, 28, 2)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 2), L[28,28],h(2, 28, 0)) ) ) ) )
Eq.ann: {}
Entry 19:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 2), L[28,28],h(2, 28, 0)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 2), L[28,28],h(1, 28, 2)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 2), L[28,28],h(2, 28, 0)) ) ) ) )
Eq.ann: {}
Entry 20:
Eq: Tile( (1, 1), G(h(1, 28, 1), L[28,28],h(1, 28, 1)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, 1), L[28,28],h(1, 28, 1)) ) )
Eq.ann: {}
Entry 21:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 2), L[28,28],h(1, 28, 1)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 2), L[28,28],h(1, 28, 1)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 1), L[28,28],h(1, 28, 1)) ) ) )
Eq.ann: {}
Entry 22:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 2), L[28,28],h(1, 28, 0)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 2), L[28,28],h(1, 28, 0)) ) ) - ( Tile( (1, 1), Tile( (4, 4), G(h(2, 28, 2), L[28,28],h(1, 28, 1)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 1), L[28,28],h(1, 28, 0)) ) ) ) )
Eq.ann: {}
Entry 23:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 1), L[28,28],h(1, 28, 0)) ) ) = Neg( ( Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 1), L[28,28],h(1, 28, 1)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 1), L[28,28],h(1, 28, 0)) ) ) ) )
Eq.ann: {}
Entry 24:
Eq: Tile( (1, 1), G(h(1, 28, 0), L[28,28],h(1, 28, 0)) ) = ( Tile( (1, 1), 1[1,1] ) Div Tile( (1, 1), G(h(1, 28, 0), L[28,28],h(1, 28, 0)) ) )
Eq.ann: {}
Entry 25:
Eq: Tile( (1, 1), Tile( (4, 4), G(h(3, 28, 1), L[28,28],h(1, 28, 0)) ) ) = ( Tile( (1, 1), Tile( (4, 4), G(h(3, 28, 1), L[28,28],h(1, 28, 0)) ) ) Kro Tile( (1, 1), Tile( (4, 4), G(h(1, 28, 0), L[28,28],h(1, 28, 0)) ) ) )
Eq.ann: {}
Entry 26:
Eq: Tile( (2, 2), Tile( (4, 4), G(h(24, 28, 4), L[28,28],h(4, 28, 0)) ) ) = ( Tile( (2, 2), Tile( (4, 4), G(h(24, 28, 4), L[28,28],h(4, 28, 0)) ) ) * Tile( (2, 2), Tile( (4, 4), G(h(4, 28, 0), L[28,28],h(4, 28, 0)) ) ) )
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
	_t0_56, _t0_57, _t0_58, _t0_59;
  __m256d _t1_0, _t1_1, _t1_2, _t1_3, _t1_4, _t1_5, _t1_6, _t1_7;
  __m256d _t2_0, _t2_1, _t2_2, _t2_3, _t2_4, _t2_5, _t2_6, _t2_7,
	_t2_8, _t2_9, _t2_10, _t2_11, _t2_12, _t2_13, _t2_14, _t2_15,
	_t2_16, _t2_17, _t2_18, _t2_19, _t2_20, _t2_21, _t2_22, _t2_23,
	_t2_24, _t2_25, _t2_26, _t2_27, _t2_28, _t2_29, _t2_30, _t2_31,
	_t2_32, _t2_33, _t2_34, _t2_35, _t2_36, _t2_37, _t2_38, _t2_39,
	_t2_40, _t2_41, _t2_42, _t2_43, _t2_44, _t2_45, _t2_46, _t2_47,
	_t2_48, _t2_49, _t2_50, _t2_51, _t2_52, _t2_53, _t2_54, _t2_55,
	_t2_56, _t2_57, _t2_58, _t2_59;
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
	_t7_56, _t7_57, _t7_58, _t7_59;
  __m256d _t8_0, _t8_1, _t8_2, _t8_3, _t8_4, _t8_5, _t8_6, _t8_7,
	_t8_8, _t8_9, _t8_10, _t8_11, _t8_12, _t8_13, _t8_14, _t8_15,
	_t8_16, _t8_17, _t8_18, _t8_19;

  _t0_1 = _mm256_maskload_pd(L + 783, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
  _t0_2 = _mm256_maskload_pd(L + 780, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, (__int64)1 << 63, 0));
  _t0_3 = _mm256_maskload_pd(L + 754, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
  _t0_6 = _mm256_maskload_pd(L + 752, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, 0, 0));
  _t0_7 = _mm256_maskload_pd(L + 725, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
  _t0_0 = _mm256_broadcast_sd(&(L[724]));
  _t0_10 = _mm256_maskload_pd(L + 724, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
  _t0_11 = _mm256_maskload_pd(L + 696, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));

  // Constant 1x1 -> 1x4
  _t0_44 = _mm256_set_pd(0, 0, 0, 1);

  // 1x1 -> 1x4
  _t0_50 = _t0_1;

  // 4-BLAC: 1x4 / 1x4
  _t0_18 = _mm256_div_pd(_t0_44, _t0_50);
  _t0_1 = _t0_18;

  // 1x1 -> 1x4
  _t0_19 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_1, _t0_1, 32), _mm256_permute2f128_pd(_t0_1, _t0_1, 32), 0);

  // 1x3 -> 1x4
  _t0_20 = _t0_2;

  // 4-BLAC: 1x4 Kro 1x4
  _t0_21 = _mm256_mul_pd(_t0_19, _t0_20);

  // 4-BLAC: -( 1x4 )
  _t0_22 = _mm256_sub_pd(_mm256_setzero_pd(), _t0_21);
  _t0_2 = _t0_22;

  // Constant 1x1 -> 1x4
  _t0_23 = _mm256_set_pd(0, 0, 0, 1);

  // 1x1 -> 1x4
  _t0_24 = _t0_3;

  // 4-BLAC: 1x4 / 1x4
  _t0_25 = _mm256_div_pd(_t0_23, _t0_24);
  _t0_3 = _t0_25;

  // 1x1 -> 1x4
  _t0_26 = _mm256_permute2f128_pd(_mm256_blend_pd(_mm256_setzero_pd(), _t0_2, 4), _mm256_blend_pd(_mm256_setzero_pd(), _t0_2, 4), 129);

  // 1x1 -> 1x4
  _t0_27 = _t0_3;

  // 4-BLAC: 1x4 Kro 1x4
  _t0_28 = _mm256_mul_pd(_t0_26, _t0_27);
  _t0_4 = _t0_28;

  // 1x2 -> 1x4
  _t0_29 = _mm256_blend_pd(_mm256_setzero_pd(), _t0_2, 3);

  // 1x1 -> 1x4
  _t0_30 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_4, _t0_4, 32), _mm256_permute2f128_pd(_t0_4, _t0_4, 32), 0);

  // 1x2 -> 1x4
  _t0_31 = _t0_6;

  // 4-BLAC: 1x4 Kro 1x4
  _t0_32 = _mm256_mul_pd(_t0_30, _t0_31);

  // 4-BLAC: 1x4 - 1x4
  _t0_33 = _mm256_sub_pd(_t0_29, _t0_32);
  _t0_5 = _t0_33;

  // 1x1 -> 1x4
  _t0_34 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_3, _t0_3, 32), _mm256_permute2f128_pd(_t0_3, _t0_3, 32), 0);

  // 1x2 -> 1x4
  _t0_35 = _t0_6;

  // 4-BLAC: 1x4 Kro 1x4
  _t0_36 = _mm256_mul_pd(_t0_34, _t0_35);

  // 4-BLAC: -( 1x4 )
  _t0_37 = _mm256_sub_pd(_mm256_setzero_pd(), _t0_36);
  _t0_6 = _t0_37;

  // Constant 1x1 -> 1x4
  _t0_38 = _mm256_set_pd(0, 0, 0, 1);

  // 1x1 -> 1x4
  _t0_39 = _t0_7;

  // 4-BLAC: 1x4 / 1x4
  _t0_40 = _mm256_div_pd(_t0_38, _t0_39);
  _t0_7 = _t0_40;

  // 2x1 -> 4x1
  _t0_41 = _mm256_unpackhi_pd(_mm256_blend_pd(_t0_6, _mm256_setzero_pd(), 12), _mm256_blend_pd(_t0_5, _mm256_setzero_pd(), 12));

  // 1x1 -> 1x4
  _t0_42 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_7, _t0_7, 32), _mm256_permute2f128_pd(_t0_7, _t0_7, 32), 0);

  // 4-BLAC: 4x1 Kro 1x4
  _t0_43 = _mm256_mul_pd(_t0_41, _t0_42);
  _t0_8 = _t0_43;

  // 2x1 -> 4x1
  _t0_45 = _mm256_blend_pd(_mm256_unpacklo_pd(_t0_6, _t0_5), _mm256_setzero_pd(), 12);

  // 2x1 -> 4x1
  _t0_46 = _t0_8;

  // 1x1 -> 1x4
  _t0_47 = _t0_0;

  // 4-BLAC: 4x1 Kro 1x4
  _t0_48 = _mm256_mul_pd(_t0_46, _t0_47);

  // 4-BLAC: 4x1 - 4x1
  _t0_49 = _mm256_sub_pd(_t0_45, _t0_48);
  _t0_9 = _t0_49;

  // 1x1 -> 1x4
  _t0_51 = _t0_7;

  // 1x1 -> 1x4
  _t0_52 = _t0_10;

  // 4-BLAC: 1x4 Kro 1x4
  _t0_53 = _mm256_mul_pd(_t0_51, _t0_52);

  // 4-BLAC: -( 1x4 )
  _t0_54 = _mm256_sub_pd(_mm256_setzero_pd(), _t0_53);
  _t0_10 = _t0_54;

  // Constant 1x1 -> 1x4
  _t0_55 = _mm256_set_pd(0, 0, 0, 1);

  // 1x1 -> 1x4
  _t0_56 = _t0_11;

  // 4-BLAC: 1x4 / 1x4
  _t0_57 = _mm256_div_pd(_t0_55, _t0_56);
  _t0_11 = _t0_57;

  // 3x1 -> 4x1
  _t0_58 = _mm256_blend_pd(_mm256_unpacklo_pd(_t0_10, _t0_9), _mm256_unpackhi_pd(_mm256_permute2f128_pd(_t0_9, _t0_9, 8), _mm256_setzero_pd()), 12);

  // 1x1 -> 1x4
  _t0_59 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_11, _t0_11, 32), _mm256_permute2f128_pd(_t0_11, _t0_11, 32), 0);

  // 4-BLAC: 4x1 Kro 1x4
  _t0_13 = _mm256_mul_pd(_t0_58, _t0_59);
  _t0_12 = _t0_13;

  // 4x4 -> 4x4 - LowTriang
  _t0_14 = _t0_11;
  _t0_15 = _mm256_blend_pd(_mm256_unpacklo_pd(_t0_12, _t0_7), _mm256_setzero_pd(), 12);
  _t0_16 = _mm256_permute2f128_pd(_mm256_shuffle_pd(_t0_12, _t0_8, 1), _mm256_unpacklo_pd(_t0_3, _mm256_setzero_pd()), 32);
  _t0_17 = _mm256_blend_pd(_mm256_permute2f128_pd(_t0_12, _mm256_unpacklo_pd(_t0_4, _t0_1), 33), _t0_8, 2);


  for( int j49 = 0; j49 <= 23; j49+=4 ) {
    _t1_4 = _mm256_loadu_pd(L + j49 + 672);
    _t1_5 = _mm256_loadu_pd(L + j49 + 700);
    _t1_6 = _mm256_loadu_pd(L + j49 + 728);
    _t1_7 = _mm256_loadu_pd(L + j49 + 756);

    // 4-BLAC: 4x4 * 4x4
    _t1_0 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_14, _t0_14, 32), _mm256_permute2f128_pd(_t0_14, _t0_14, 32), 0), _t1_4), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_14, _t0_14, 32), _mm256_permute2f128_pd(_t0_14, _t0_14, 32), 15), _t1_5)), _mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_14, _t0_14, 49), _mm256_permute2f128_pd(_t0_14, _t0_14, 49), 0), _t1_6), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_14, _t0_14, 49), _mm256_permute2f128_pd(_t0_14, _t0_14, 49), 15), _t1_7)));
    _t1_1 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_15, _t0_15, 32), _mm256_permute2f128_pd(_t0_15, _t0_15, 32), 0), _t1_4), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_15, _t0_15, 32), _mm256_permute2f128_pd(_t0_15, _t0_15, 32), 15), _t1_5)), _mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_15, _t0_15, 49), _mm256_permute2f128_pd(_t0_15, _t0_15, 49), 0), _t1_6), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_15, _t0_15, 49), _mm256_permute2f128_pd(_t0_15, _t0_15, 49), 15), _t1_7)));
    _t1_2 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_16, _t0_16, 32), _mm256_permute2f128_pd(_t0_16, _t0_16, 32), 0), _t1_4), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_16, _t0_16, 32), _mm256_permute2f128_pd(_t0_16, _t0_16, 32), 15), _t1_5)), _mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_16, _t0_16, 49), _mm256_permute2f128_pd(_t0_16, _t0_16, 49), 0), _t1_6), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_16, _t0_16, 49), _mm256_permute2f128_pd(_t0_16, _t0_16, 49), 15), _t1_7)));
    _t1_3 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_17, _t0_17, 32), _mm256_permute2f128_pd(_t0_17, _t0_17, 32), 0), _t1_4), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_17, _t0_17, 32), _mm256_permute2f128_pd(_t0_17, _t0_17, 32), 15), _t1_5)), _mm256_add_pd(_mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_17, _t0_17, 49), _mm256_permute2f128_pd(_t0_17, _t0_17, 49), 0), _t1_6), _mm256_mul_pd(_mm256_shuffle_pd(_mm256_permute2f128_pd(_t0_17, _t0_17, 49), _mm256_permute2f128_pd(_t0_17, _t0_17, 49), 15), _t1_7)));

    // 4-BLAC: -( 4x4 )
    _t1_4 = _mm256_sub_pd(_mm256_setzero_pd(), _t1_0);
    _t1_5 = _mm256_sub_pd(_mm256_setzero_pd(), _t1_1);
    _t1_6 = _mm256_sub_pd(_mm256_setzero_pd(), _t1_2);
    _t1_7 = _mm256_sub_pd(_mm256_setzero_pd(), _t1_3);
    _mm256_storeu_pd(L + j49 + 672, _t1_4);
    _mm256_storeu_pd(L + j49 + 700, _t1_5);
    _mm256_storeu_pd(L + j49 + 728, _t1_6);
    _mm256_storeu_pd(L + j49 + 756, _t1_7);
  }


  for( int fi2 = 4; fi2 <= 23; fi2+=4 ) {
    _t2_1 = _mm256_maskload_pd(L + -29*fi2 + 783, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
    _t2_2 = _mm256_maskload_pd(L + -29*fi2 + 780, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, (__int64)1 << 63, 0));
    _t2_3 = _mm256_maskload_pd(L + -29*fi2 + 754, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
    _t2_6 = _mm256_maskload_pd(L + -29*fi2 + 752, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, 0, 0));
    _t2_7 = _mm256_maskload_pd(L + -29*fi2 + 725, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
    _t2_0 = _mm256_broadcast_sd(&(L[-29*fi2 + 724]));
    _t2_10 = _mm256_maskload_pd(L + -29*fi2 + 724, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
    _t2_11 = _mm256_maskload_pd(L + -29*fi2 + 696, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));

    // Constant 1x1 -> 1x4
    _t2_13 = _mm256_set_pd(0, 0, 0, 1);

    // 1x1 -> 1x4
    _t2_14 = _t2_1;

    // 4-BLAC: 1x4 / 1x4
    _t2_15 = _mm256_div_pd(_t2_13, _t2_14);
    _t2_1 = _t2_15;

    // 1x1 -> 1x4
    _t2_16 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t2_1, _t2_1, 32), _mm256_permute2f128_pd(_t2_1, _t2_1, 32), 0);

    // 1x3 -> 1x4
    _t2_17 = _t2_2;

    // 4-BLAC: 1x4 Kro 1x4
    _t2_18 = _mm256_mul_pd(_t2_16, _t2_17);

    // 4-BLAC: -( 1x4 )
    _t2_19 = _mm256_sub_pd(_mm256_setzero_pd(), _t2_18);
    _t2_2 = _t2_19;

    // Constant 1x1 -> 1x4
    _t2_20 = _mm256_set_pd(0, 0, 0, 1);

    // 1x1 -> 1x4
    _t2_21 = _t2_3;

    // 4-BLAC: 1x4 / 1x4
    _t2_22 = _mm256_div_pd(_t2_20, _t2_21);
    _t2_3 = _t2_22;

    // 1x1 -> 1x4
    _t2_23 = _mm256_permute2f128_pd(_mm256_blend_pd(_mm256_setzero_pd(), _t2_2, 4), _mm256_blend_pd(_mm256_setzero_pd(), _t2_2, 4), 129);

    // 1x1 -> 1x4
    _t2_24 = _t2_3;

    // 4-BLAC: 1x4 Kro 1x4
    _t2_25 = _mm256_mul_pd(_t2_23, _t2_24);
    _t2_4 = _t2_25;

    // 1x2 -> 1x4
    _t2_26 = _mm256_blend_pd(_mm256_setzero_pd(), _t2_2, 3);

    // 1x1 -> 1x4
    _t2_27 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t2_4, _t2_4, 32), _mm256_permute2f128_pd(_t2_4, _t2_4, 32), 0);

    // 1x2 -> 1x4
    _t2_28 = _t2_6;

    // 4-BLAC: 1x4 Kro 1x4
    _t2_29 = _mm256_mul_pd(_t2_27, _t2_28);

    // 4-BLAC: 1x4 - 1x4
    _t2_30 = _mm256_sub_pd(_t2_26, _t2_29);
    _t2_5 = _t2_30;

    // 1x1 -> 1x4
    _t2_31 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t2_3, _t2_3, 32), _mm256_permute2f128_pd(_t2_3, _t2_3, 32), 0);

    // 1x2 -> 1x4
    _t2_32 = _t2_6;

    // 4-BLAC: 1x4 Kro 1x4
    _t2_33 = _mm256_mul_pd(_t2_31, _t2_32);

    // 4-BLAC: -( 1x4 )
    _t2_34 = _mm256_sub_pd(_mm256_setzero_pd(), _t2_33);
    _t2_6 = _t2_34;

    // Constant 1x1 -> 1x4
    _t2_35 = _mm256_set_pd(0, 0, 0, 1);

    // 1x1 -> 1x4
    _t2_36 = _t2_7;

    // 4-BLAC: 1x4 / 1x4
    _t2_37 = _mm256_div_pd(_t2_35, _t2_36);
    _t2_7 = _t2_37;

    // 2x1 -> 4x1
    _t2_38 = _mm256_unpackhi_pd(_mm256_blend_pd(_t2_6, _mm256_setzero_pd(), 12), _mm256_blend_pd(_t2_5, _mm256_setzero_pd(), 12));

    // 1x1 -> 1x4
    _t2_39 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t2_7, _t2_7, 32), _mm256_permute2f128_pd(_t2_7, _t2_7, 32), 0);

    // 4-BLAC: 4x1 Kro 1x4
    _t2_40 = _mm256_mul_pd(_t2_38, _t2_39);
    _t2_8 = _t2_40;

    // 2x1 -> 4x1
    _t2_41 = _mm256_blend_pd(_mm256_unpacklo_pd(_t2_6, _t2_5), _mm256_setzero_pd(), 12);

    // 2x1 -> 4x1
    _t2_42 = _t2_8;

    // 1x1 -> 1x4
    _t2_43 = _t2_0;

    // 4-BLAC: 4x1 Kro 1x4
    _t2_44 = _mm256_mul_pd(_t2_42, _t2_43);

    // 4-BLAC: 4x1 - 4x1
    _t2_45 = _mm256_sub_pd(_t2_41, _t2_44);
    _t2_9 = _t2_45;

    // 1x1 -> 1x4
    _t2_46 = _t2_7;

    // 1x1 -> 1x4
    _t2_47 = _t2_10;

    // 4-BLAC: 1x4 Kro 1x4
    _t2_48 = _mm256_mul_pd(_t2_46, _t2_47);

    // 4-BLAC: -( 1x4 )
    _t2_49 = _mm256_sub_pd(_mm256_setzero_pd(), _t2_48);
    _t2_10 = _t2_49;

    // Constant 1x1 -> 1x4
    _t2_50 = _mm256_set_pd(0, 0, 0, 1);

    // 1x1 -> 1x4
    _t2_51 = _t2_11;

    // 4-BLAC: 1x4 / 1x4
    _t2_52 = _mm256_div_pd(_t2_50, _t2_51);
    _t2_11 = _t2_52;

    // 3x1 -> 4x1
    _t2_53 = _mm256_blend_pd(_mm256_unpacklo_pd(_t2_10, _t2_9), _mm256_unpackhi_pd(_mm256_permute2f128_pd(_t2_9, _t2_9, 8), _mm256_setzero_pd()), 12);

    // 1x1 -> 1x4
    _t2_54 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t2_11, _t2_11, 32), _mm256_permute2f128_pd(_t2_11, _t2_11, 32), 0);

    // 4-BLAC: 4x1 Kro 1x4
    _t2_55 = _mm256_mul_pd(_t2_53, _t2_54);
    _t2_12 = _t2_55;

    // 4x4 -> 4x4 - LowTriang
    _t2_56 = _t2_11;
    _t2_57 = _mm256_blend_pd(_mm256_unpacklo_pd(_t2_12, _t2_7), _mm256_setzero_pd(), 12);
    _t2_58 = _mm256_permute2f128_pd(_mm256_shuffle_pd(_t2_12, _t2_8, 1), _mm256_unpacklo_pd(_t2_3, _mm256_setzero_pd()), 32);
    _t2_59 = _mm256_blend_pd(_mm256_permute2f128_pd(_t2_12, _mm256_unpacklo_pd(_t2_4, _t2_1), 33), _t2_8, 2);
    _mm256_maskstore_pd(L + -29*fi2 + 780, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, (__int64)1 << 63, 0), _t2_2);
    _mm256_maskstore_pd(L + -29*fi2 + 780, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, 0, 0), _t2_5);
    _mm256_maskstore_pd(L + -29*fi2 + 752, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, 0, 0), _t2_6);
    _mm256_maskstore_pd(L + -29*fi2 + 752, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_9);
    _mm256_maskstore_pd(L + -29*fi2 + 780, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t2_9, _t2_9, 1));
    _mm256_maskstore_pd(L + -29*fi2 + 724, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_10);

    for( int j49 = 0; j49 <= fi2 - 1; j49+=4 ) {
      _t3_15 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 808);
      _t3_14 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 809);
      _t3_13 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 810);
      _t3_12 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 811);
      _t3_11 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 836);
      _t3_10 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 837);
      _t3_9 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 838);
      _t3_8 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 839);
      _t3_7 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 864);
      _t3_6 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 865);
      _t3_5 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 866);
      _t3_4 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 867);
      _t3_3 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 892);
      _t3_2 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 893);
      _t3_1 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 894);
      _t3_0 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 895);

      // 4x4 -> 4x4 - LowTriang
      _t3_20 = _t2_11;
      _t3_21 = _mm256_blend_pd(_mm256_unpacklo_pd(_t2_12, _t2_7), _mm256_setzero_pd(), 12);
      _t3_22 = _mm256_permute2f128_pd(_mm256_shuffle_pd(_t2_12, _t2_8, 1), _mm256_unpacklo_pd(_t2_3, _mm256_setzero_pd()), 32);
      _t3_23 = _mm256_blend_pd(_mm256_permute2f128_pd(_t2_12, _mm256_unpacklo_pd(_t2_4, _t2_1), 33), _t2_8, 2);

      // 4-BLAC: 4x4 * 4x4
      _t3_16 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t3_15, _t3_20), _mm256_mul_pd(_t3_14, _t3_21)), _mm256_add_pd(_mm256_mul_pd(_t3_13, _t3_22), _mm256_mul_pd(_t3_12, _t3_23)));
      _t3_17 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t3_11, _t3_20), _mm256_mul_pd(_t3_10, _t3_21)), _mm256_add_pd(_mm256_mul_pd(_t3_9, _t3_22), _mm256_mul_pd(_t3_8, _t3_23)));
      _t3_18 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t3_7, _t3_20), _mm256_mul_pd(_t3_6, _t3_21)), _mm256_add_pd(_mm256_mul_pd(_t3_5, _t3_22), _mm256_mul_pd(_t3_4, _t3_23)));
      _t3_19 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t3_3, _t3_20), _mm256_mul_pd(_t3_2, _t3_21)), _mm256_add_pd(_mm256_mul_pd(_t3_1, _t3_22), _mm256_mul_pd(_t3_0, _t3_23)));
      _mm256_storeu_pd(L + -29*fi2 + 28*j49 + 808, _t3_16);
      _mm256_storeu_pd(L + -29*fi2 + 28*j49 + 836, _t3_17);
      _mm256_storeu_pd(L + -29*fi2 + 28*j49 + 864, _t3_18);
      _mm256_storeu_pd(L + -29*fi2 + 28*j49 + 892, _t3_19);
    }

    for( int j49 = 0; j49 <= fi2 - 1; j49+=4 ) {

      for( int i99 = 0; i99 <= -fi2 + 23; i99+=4 ) {
        _t4_24 = _mm256_loadu_pd(L + -28*fi2 + i99 + 28*j49 + 784);
        _t4_25 = _mm256_loadu_pd(L + -28*fi2 + i99 + 28*j49 + 812);
        _t4_26 = _mm256_loadu_pd(L + -28*fi2 + i99 + 28*j49 + 840);
        _t4_27 = _mm256_loadu_pd(L + -28*fi2 + i99 + 28*j49 + 868);
        _t4_19 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 808);
        _t4_18 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 809);
        _t4_17 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 810);
        _t4_16 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 811);
        _t4_15 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 836);
        _t4_14 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 837);
        _t4_13 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 838);
        _t4_12 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 839);
        _t4_11 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 864);
        _t4_10 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 865);
        _t4_9 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 866);
        _t4_8 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 867);
        _t4_7 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 892);
        _t4_6 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 893);
        _t4_5 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 894);
        _t4_4 = _mm256_broadcast_sd(L + -29*fi2 + 28*j49 + 895);
        _t4_3 = _mm256_loadu_pd(L + -28*fi2 + i99 + 672);
        _t4_2 = _mm256_loadu_pd(L + -28*fi2 + i99 + 700);
        _t4_1 = _mm256_loadu_pd(L + -28*fi2 + i99 + 728);
        _t4_0 = _mm256_loadu_pd(L + -28*fi2 + i99 + 756);

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
        _mm256_storeu_pd(L + -28*fi2 + i99 + 28*j49 + 784, _t4_24);
        _mm256_storeu_pd(L + -28*fi2 + i99 + 28*j49 + 812, _t4_25);
        _mm256_storeu_pd(L + -28*fi2 + i99 + 28*j49 + 840, _t4_26);
        _mm256_storeu_pd(L + -28*fi2 + i99 + 28*j49 + 868, _t4_27);
      }
    }

    // 4x4 -> 4x4 - LowTriang
    _t5_0 = _t2_11;
    _t5_1 = _mm256_blend_pd(_mm256_unpacklo_pd(_t2_12, _t2_7), _mm256_setzero_pd(), 12);
    _t5_2 = _mm256_permute2f128_pd(_mm256_shuffle_pd(_t2_12, _t2_8, 1), _mm256_unpacklo_pd(_t2_3, _mm256_setzero_pd()), 32);
    _t5_3 = _mm256_blend_pd(_mm256_permute2f128_pd(_t2_12, _mm256_unpacklo_pd(_t2_4, _t2_1), 33), _t2_8, 2);

    for( int j49 = 0; j49 <= -fi2 + 23; j49+=4 ) {
      _t6_4 = _mm256_loadu_pd(L + -28*fi2 + j49 + 672);
      _t6_5 = _mm256_loadu_pd(L + -28*fi2 + j49 + 700);
      _t6_6 = _mm256_loadu_pd(L + -28*fi2 + j49 + 728);
      _t6_7 = _mm256_loadu_pd(L + -28*fi2 + j49 + 756);

      // 4x4 -> 4x4 - LowTriang
      _t6_8 = _t2_11;
      _t6_9 = _mm256_blend_pd(_mm256_unpacklo_pd(_t2_12, _t2_7), _mm256_setzero_pd(), 12);
      _t6_10 = _mm256_permute2f128_pd(_mm256_shuffle_pd(_t2_12, _t2_8, 1), _mm256_unpacklo_pd(_t2_3, _mm256_setzero_pd()), 32);
      _t6_11 = _mm256_blend_pd(_mm256_permute2f128_pd(_t2_12, _mm256_unpacklo_pd(_t2_4, _t2_1), 33), _t2_8, 2);

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
      _mm256_storeu_pd(L + -28*fi2 + j49 + 672, _t6_4);
      _mm256_storeu_pd(L + -28*fi2 + j49 + 700, _t6_5);
      _mm256_storeu_pd(L + -28*fi2 + j49 + 728, _t6_6);
      _mm256_storeu_pd(L + -28*fi2 + j49 + 756, _t6_7);
    }
    _mm256_maskstore_pd(L + -29*fi2 + 783, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_1);
    _mm256_maskstore_pd(L + -29*fi2 + 754, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_3);
    _mm256_maskstore_pd(L + -29*fi2 + 782, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_4);
    _mm256_maskstore_pd(L + -29*fi2 + 725, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_7);
    _mm256_maskstore_pd(L + -29*fi2 + 753, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_8);
    _mm256_maskstore_pd(L + -29*fi2 + 781, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t2_8, _t2_8, 1));
    _mm256_maskstore_pd(L + -29*fi2 + 696, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_11);
    _mm256_maskstore_pd(L + -29*fi2 + 724, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t2_12);
    _mm256_maskstore_pd(L + -29*fi2 + 752, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t2_12, _t2_12, 1));
    _mm256_maskstore_pd(L + -29*fi2 + 780, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_permute2f128_pd(_t2_12, _t2_12, 129));
  }

  _t7_1 = _mm256_maskload_pd(L + 87, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
  _t7_2 = _mm256_maskload_pd(L + 84, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, (__int64)1 << 63, 0));
  _t7_3 = _mm256_maskload_pd(L + 58, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
  _t7_6 = _mm256_maskload_pd(L + 56, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, 0, 0));
  _t7_7 = _mm256_maskload_pd(L + 29, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
  _t7_0 = _mm256_broadcast_sd(&(L[28]));
  _t7_10 = _mm256_maskload_pd(L + 28, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));
  _t7_11 = _mm256_maskload_pd(L, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0));

  // Constant 1x1 -> 1x4
  _t7_13 = _mm256_set_pd(0, 0, 0, 1);

  // 1x1 -> 1x4
  _t7_14 = _t7_1;

  // 4-BLAC: 1x4 / 1x4
  _t7_15 = _mm256_div_pd(_t7_13, _t7_14);
  _t7_1 = _t7_15;

  // 1x1 -> 1x4
  _t7_16 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_1, _t7_1, 32), _mm256_permute2f128_pd(_t7_1, _t7_1, 32), 0);

  // 1x3 -> 1x4
  _t7_17 = _t7_2;

  // 4-BLAC: 1x4 Kro 1x4
  _t7_18 = _mm256_mul_pd(_t7_16, _t7_17);

  // 4-BLAC: -( 1x4 )
  _t7_19 = _mm256_sub_pd(_mm256_setzero_pd(), _t7_18);
  _t7_2 = _t7_19;

  // Constant 1x1 -> 1x4
  _t7_20 = _mm256_set_pd(0, 0, 0, 1);

  // 1x1 -> 1x4
  _t7_21 = _t7_3;

  // 4-BLAC: 1x4 / 1x4
  _t7_22 = _mm256_div_pd(_t7_20, _t7_21);
  _t7_3 = _t7_22;

  // 1x1 -> 1x4
  _t7_23 = _mm256_permute2f128_pd(_mm256_blend_pd(_mm256_setzero_pd(), _t7_2, 4), _mm256_blend_pd(_mm256_setzero_pd(), _t7_2, 4), 129);

  // 1x1 -> 1x4
  _t7_24 = _t7_3;

  // 4-BLAC: 1x4 Kro 1x4
  _t7_25 = _mm256_mul_pd(_t7_23, _t7_24);
  _t7_4 = _t7_25;

  // 1x2 -> 1x4
  _t7_26 = _mm256_blend_pd(_mm256_setzero_pd(), _t7_2, 3);

  // 1x1 -> 1x4
  _t7_27 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_4, _t7_4, 32), _mm256_permute2f128_pd(_t7_4, _t7_4, 32), 0);

  // 1x2 -> 1x4
  _t7_28 = _t7_6;

  // 4-BLAC: 1x4 Kro 1x4
  _t7_29 = _mm256_mul_pd(_t7_27, _t7_28);

  // 4-BLAC: 1x4 - 1x4
  _t7_30 = _mm256_sub_pd(_t7_26, _t7_29);
  _t7_5 = _t7_30;

  // 1x1 -> 1x4
  _t7_31 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_3, _t7_3, 32), _mm256_permute2f128_pd(_t7_3, _t7_3, 32), 0);

  // 1x2 -> 1x4
  _t7_32 = _t7_6;

  // 4-BLAC: 1x4 Kro 1x4
  _t7_33 = _mm256_mul_pd(_t7_31, _t7_32);

  // 4-BLAC: -( 1x4 )
  _t7_34 = _mm256_sub_pd(_mm256_setzero_pd(), _t7_33);
  _t7_6 = _t7_34;

  // Constant 1x1 -> 1x4
  _t7_35 = _mm256_set_pd(0, 0, 0, 1);

  // 1x1 -> 1x4
  _t7_36 = _t7_7;

  // 4-BLAC: 1x4 / 1x4
  _t7_37 = _mm256_div_pd(_t7_35, _t7_36);
  _t7_7 = _t7_37;

  // 2x1 -> 4x1
  _t7_38 = _mm256_unpackhi_pd(_mm256_blend_pd(_t7_6, _mm256_setzero_pd(), 12), _mm256_blend_pd(_t7_5, _mm256_setzero_pd(), 12));

  // 1x1 -> 1x4
  _t7_39 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_7, _t7_7, 32), _mm256_permute2f128_pd(_t7_7, _t7_7, 32), 0);

  // 4-BLAC: 4x1 Kro 1x4
  _t7_40 = _mm256_mul_pd(_t7_38, _t7_39);
  _t7_8 = _t7_40;

  // 2x1 -> 4x1
  _t7_41 = _mm256_blend_pd(_mm256_unpacklo_pd(_t7_6, _t7_5), _mm256_setzero_pd(), 12);

  // 2x1 -> 4x1
  _t7_42 = _t7_8;

  // 1x1 -> 1x4
  _t7_43 = _t7_0;

  // 4-BLAC: 4x1 Kro 1x4
  _t7_44 = _mm256_mul_pd(_t7_42, _t7_43);

  // 4-BLAC: 4x1 - 4x1
  _t7_45 = _mm256_sub_pd(_t7_41, _t7_44);
  _t7_9 = _t7_45;

  // 1x1 -> 1x4
  _t7_46 = _t7_7;

  // 1x1 -> 1x4
  _t7_47 = _t7_10;

  // 4-BLAC: 1x4 Kro 1x4
  _t7_48 = _mm256_mul_pd(_t7_46, _t7_47);

  // 4-BLAC: -( 1x4 )
  _t7_49 = _mm256_sub_pd(_mm256_setzero_pd(), _t7_48);
  _t7_10 = _t7_49;

  // Constant 1x1 -> 1x4
  _t7_50 = _mm256_set_pd(0, 0, 0, 1);

  // 1x1 -> 1x4
  _t7_51 = _t7_11;

  // 4-BLAC: 1x4 / 1x4
  _t7_52 = _mm256_div_pd(_t7_50, _t7_51);
  _t7_11 = _t7_52;

  // 3x1 -> 4x1
  _t7_53 = _mm256_blend_pd(_mm256_unpacklo_pd(_t7_10, _t7_9), _mm256_unpackhi_pd(_mm256_permute2f128_pd(_t7_9, _t7_9, 8), _mm256_setzero_pd()), 12);

  // 1x1 -> 1x4
  _t7_54 = _mm256_shuffle_pd(_mm256_permute2f128_pd(_t7_11, _t7_11, 32), _mm256_permute2f128_pd(_t7_11, _t7_11, 32), 0);

  // 4-BLAC: 4x1 Kro 1x4
  _t7_55 = _mm256_mul_pd(_t7_53, _t7_54);
  _t7_12 = _t7_55;

  // 4x4 -> 4x4 - LowTriang
  _t7_56 = _t7_11;
  _t7_57 = _mm256_blend_pd(_mm256_unpacklo_pd(_t7_12, _t7_7), _mm256_setzero_pd(), 12);
  _t7_58 = _mm256_permute2f128_pd(_mm256_shuffle_pd(_t7_12, _t7_8, 1), _mm256_unpacklo_pd(_t7_3, _mm256_setzero_pd()), 32);
  _t7_59 = _mm256_blend_pd(_mm256_permute2f128_pd(_t7_12, _mm256_unpacklo_pd(_t7_4, _t7_1), 33), _t7_8, 2);


  for( int j49 = 0; j49 <= 23; j49+=8 ) {

    for( int i99 = 0; i99 <= 7; i99+=4 ) {
      _t8_15 = _mm256_broadcast_sd(L + 28*i99 + 28*j49 + 112);
      _t8_14 = _mm256_broadcast_sd(L + 28*i99 + 28*j49 + 113);
      _t8_13 = _mm256_broadcast_sd(L + 28*i99 + 28*j49 + 114);
      _t8_12 = _mm256_broadcast_sd(L + 28*i99 + 28*j49 + 115);
      _t8_11 = _mm256_broadcast_sd(L + 28*i99 + 28*j49 + 140);
      _t8_10 = _mm256_broadcast_sd(L + 28*i99 + 28*j49 + 141);
      _t8_9 = _mm256_broadcast_sd(L + 28*i99 + 28*j49 + 142);
      _t8_8 = _mm256_broadcast_sd(L + 28*i99 + 28*j49 + 143);
      _t8_7 = _mm256_broadcast_sd(L + 28*i99 + 28*j49 + 168);
      _t8_6 = _mm256_broadcast_sd(L + 28*i99 + 28*j49 + 169);
      _t8_5 = _mm256_broadcast_sd(L + 28*i99 + 28*j49 + 170);
      _t8_4 = _mm256_broadcast_sd(L + 28*i99 + 28*j49 + 171);
      _t8_3 = _mm256_broadcast_sd(L + 28*i99 + 28*j49 + 196);
      _t8_2 = _mm256_broadcast_sd(L + 28*i99 + 28*j49 + 197);
      _t8_1 = _mm256_broadcast_sd(L + 28*i99 + 28*j49 + 198);
      _t8_0 = _mm256_broadcast_sd(L + 28*i99 + 28*j49 + 199);

      // 4-BLAC: 4x4 * 4x4
      _t8_16 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t8_15, _t7_56), _mm256_mul_pd(_t8_14, _t7_57)), _mm256_add_pd(_mm256_mul_pd(_t8_13, _t7_58), _mm256_mul_pd(_t8_12, _t7_59)));
      _t8_17 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t8_11, _t7_56), _mm256_mul_pd(_t8_10, _t7_57)), _mm256_add_pd(_mm256_mul_pd(_t8_9, _t7_58), _mm256_mul_pd(_t8_8, _t7_59)));
      _t8_18 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t8_7, _t7_56), _mm256_mul_pd(_t8_6, _t7_57)), _mm256_add_pd(_mm256_mul_pd(_t8_5, _t7_58), _mm256_mul_pd(_t8_4, _t7_59)));
      _t8_19 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_t8_3, _t7_56), _mm256_mul_pd(_t8_2, _t7_57)), _mm256_add_pd(_mm256_mul_pd(_t8_1, _t7_58), _mm256_mul_pd(_t8_0, _t7_59)));
      _mm256_storeu_pd(L + 28*i99 + 28*j49 + 112, _t8_16);
      _mm256_storeu_pd(L + 28*i99 + 28*j49 + 140, _t8_17);
      _mm256_storeu_pd(L + 28*i99 + 28*j49 + 168, _t8_18);
      _mm256_storeu_pd(L + 28*i99 + 28*j49 + 196, _t8_19);
    }
  }

  _mm256_maskstore_pd(L + 783, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_1);
  _mm256_maskstore_pd(L + 780, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, (__int64)1 << 63, 0), _t0_2);
  _mm256_maskstore_pd(L + 754, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_3);
  _mm256_maskstore_pd(L + 782, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_4);
  _mm256_maskstore_pd(L + 780, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, 0, 0), _t0_5);
  _mm256_maskstore_pd(L + 752, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, 0, 0), _t0_6);
  _mm256_maskstore_pd(L + 725, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_7);
  _mm256_maskstore_pd(L + 753, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_8);
  _mm256_maskstore_pd(L + 781, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t0_8, _t0_8, 1));
  _mm256_maskstore_pd(L + 752, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_9);
  _mm256_maskstore_pd(L + 780, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t0_9, _t0_9, 1));
  _mm256_maskstore_pd(L + 724, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_10);
  _mm256_maskstore_pd(L + 696, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_11);
  _mm256_maskstore_pd(L + 724, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t0_12);
  _mm256_maskstore_pd(L + 752, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t0_12, _t0_12, 1));
  _mm256_maskstore_pd(L + 780, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_permute2f128_pd(_t0_12, _t0_12, 129));
  _mm256_maskstore_pd(L + 87, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_1);
  _mm256_maskstore_pd(L + 84, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, (__int64)1 << 63, 0), _t7_2);
  _mm256_maskstore_pd(L + 58, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_3);
  _mm256_maskstore_pd(L + 86, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_4);
  _mm256_maskstore_pd(L + 84, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, 0, 0), _t7_5);
  _mm256_maskstore_pd(L + 56, _mm256_setr_epi64x((__int64)1 << 63, (__int64)1 << 63, 0, 0), _t7_6);
  _mm256_maskstore_pd(L + 29, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_7);
  _mm256_maskstore_pd(L + 57, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_8);
  _mm256_maskstore_pd(L + 85, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t7_8, _t7_8, 1));
  _mm256_maskstore_pd(L + 56, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_9);
  _mm256_maskstore_pd(L + 84, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t7_9, _t7_9, 1));
  _mm256_maskstore_pd(L + 28, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_10);
  _mm256_maskstore_pd(L, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_11);
  _mm256_maskstore_pd(L + 28, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _t7_12);
  _mm256_maskstore_pd(L + 56, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_shuffle_pd(_t7_12, _t7_12, 1));
  _mm256_maskstore_pd(L + 84, _mm256_setr_epi64x((__int64)1 << 63, 0, 0, 0), _mm256_permute2f128_pd(_t7_12, _t7_12, 129));

}
