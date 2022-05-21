#ifndef DIAG_C
#define DIAG_C

#include <complex.h>

#ifdef QUAD
#define RealType long double
#define Re creall
#define Im cimagl
#define Conjugate conjl
#define Sqrt csqrtl
#define Abs cabsl
#define absr fabsl
#define sign copysignl
#else
#define RealType double
#define Re creal
#define Im cimag
#define Conjugate conj
#define Sqrt csqrt
#define Abs cabs
#define absr fabs
#define sign copysign
#endif

typedef const int cint;
typedef const RealType cRealType;
typedef RealType complex ComplexType;
typedef const ComplexType cComplexType;

#define Element(M,i,j) M[(i)*ld##M+(j)]
#define A(i,j) Element(A,i,j)
#define U(i,j) Element(U,i,j)
#define V_(i,j) Element(V_,i,j)
#define W_(i,j) Element(W_,i,j)
#define A_(i,j) Element(A_,i,j)

#if UCOLS
#define UL(i,j) U(j,i)
#define VL(i,j) Element(V,j,i)
#define WL(i,j) Element(W,j,i)
#else
#define UL(i,j) U(i,j)
#define VL(i,j) Element(V,i,j)
#define WL(i,j) Element(W,i,j)
#endif

void TakagiFactor(cint n, ComplexType *A, cint ldA, RealType *d, ComplexType *U, cint ldU, cint sort);
  
  /* A matrix is considered diagonal if the sum of the squares
   of the off-diagonal elements is less than EPS.  SYM_EPS is
   half of EPS since only the upper triangle is counted for
   symmetric matrices.
   52 bits is the mantissa length for IEEE double precision. */

#define EPS 0x1p-102
#define SYM_EPS 0x1p-103
#define DBL_EPS 0x1p-52

#endif
