/*
   diag-c.h
   global declarations for the Diag routines in C
   this file is part of Diag
   last modified 20 Aug 15 th
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tagaki_factor.h"

static inline RealType Sq(cComplexType x) { return Re(x*Conjugate(x)); }
static inline RealType sq(cRealType x) { return x*x; }
static inline RealType min(cRealType a, cRealType b) { return (a < b) ? a : b; }
static inline RealType max(cRealType a, cRealType b) { return (a > b) ? a : b; }
static inline int imin(cint a, cint b) { return (a < b) ? a : b; }
static inline int imax(cint a, cint b) { return (a > b) ? a : b; }

int nsweeps;

/*
   TakagiFactor.c
   computes the Takagi factorization of a complex symmetric matrix
   code adapted from the "Handbook" routines
   (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
   this file is part of the Diag library
   last modified 20 Aug 15 th
*/

/*
   TakagiFactor factorizes a complex symmetric n-by-n matrix
   Input:	n, A = n-by-n matrix, complex symmetric
		(only the upper triangle of A needs to be filled),
   Output:	d = vector of diagonal values,
		U = transformation matrix, unitary (U^-1 = U^+),
   these fulfill
	d = U^* A U^+,  A = U^T d U,  U^* A = d U  (UCOLS=0),
	d = U^+ A U^*,  A = U d U^T,  A U^* = U d  (UCOLS=1).
*/

void TakagiFactor(cint n, ComplexType *A, cint ldA,
  RealType *d, ComplexType *U, cint ldU, cint sort)
{
  int p, q;
  cRealType red = .04/(n*n*n*n);
  ComplexType ev[n][2];

  for( p = 0; p < n; ++p ) {
    ev[p][0] = 0;
    ev[p][1] = A(p,p);
  }

  for( p = 0; p < n; ++p ) {
    memset(&U(p,0), 0, n*sizeof(ComplexType));
    U(p,p) = 1;
  }

  for( nsweeps = 1; nsweeps <= 50; ++nsweeps ) {
    RealType thresh = 0;
    for( q = 1; q < n; ++q )
      for( p = 0; p < q; ++p )
        thresh += Sq(A(p,q));
    if( !(thresh > SYM_EPS) ) goto done;

    thresh = (nsweeps < 4) ? thresh*red : 0;

    for( q = 1; q < n; ++q )
      for( p = 0; p < q; ++p ) {
        cComplexType Apq = A(p,q);
        cRealType off = Sq(Apq);
        cRealType sqp = Sq(ev[p][1]);
        cRealType sqq = Sq(ev[q][1]);
        if( nsweeps > 4 && off < SYM_EPS*(sqp + sqq) )
          A(p,q) = 0;
        else if( off > thresh ) {
          RealType t, invc;
          ComplexType f;
          int j;

          t = .5*absr(sqp - sqq);
          if( t > DBL_EPS )
            f = sign(1, sqp - sqq)*
              (ev[q][1]*Conjugate(Apq) + Conjugate(ev[p][1])*Apq);
          else
            f = (sqp == 0) ? 1 : Sqrt(ev[q][1]/ev[p][1]);
          t += sqrt(t*t + Sq(f));
          f /= t;

          ev[p][1] = A(p,p) + (ev[p][0] += Apq*Conjugate(f));
          ev[q][1] = A(q,q) + (ev[q][0] -= Apq*f);

          t = Sq(f);
          invc = sqrt(t + 1);
          f /= invc;
          t /= invc*(invc + 1);

          for( j = 0; j < p; ++j ) {
            cComplexType x = A(j,p);
            cComplexType y = A(j,q);
            A(j,p) = x + (Conjugate(f)*y - t*x);
            A(j,q) = y - (f*x + t*y);
          }

          for( j = p + 1; j < q; ++j ) {
            cComplexType x = A(p,j);
            cComplexType y = A(j,q);
            A(p,j) = x + (Conjugate(f)*y - t*x);
            A(j,q) = y - (f*x + t*y);
          }

          for( j = q + 1; j < n; ++j ) {
            cComplexType x = A(p,j);
            cComplexType y = A(q,j);
            A(p,j) = x + (Conjugate(f)*y - t*x);
            A(q,j) = y - (f*x + t*y);
          }

          A(p,q) = 0;

          for( j = 0; j < n; ++j ) {
            cComplexType x = UL(p,j);
            cComplexType y = UL(q,j);
            UL(p,j) = x + (f*y - t*x);
            UL(q,j) = y - (Conjugate(f)*x + t*y);
          }
        }
      }

    for( p = 0; p < n; ++p ) {
      ev[p][0] = 0;
      A(p,p) = ev[p][1];
    }
  }

  fputs("Bad convergence in TakagiFactor\n", stderr);

done:

/* make the diagonal elements nonnegative */

  for( p = 0; p < n; ++p ) {
    cComplexType App = A(p,p);
    d[p] = Abs(App);
    if( d[p] > DBL_EPS && d[p] != Re(App) ) {
      cComplexType f = Sqrt(App/d[p]);
      for( q = 0; q < n; ++q ) UL(p,q) *= f;
    }
  }

  if( sort == 0 ) return;

/* sort the eigenvalues */

  for( p = 0; p < n - 1; ++p ) {
    int j = p;
    RealType t = d[p];
    for( q = p + 1; q < n; ++q )
      if( sort*(t - d[q]) > 0 ) t = d[j = q];
    if( j == p ) continue;
    d[j] = d[p];
    d[p] = t;
    for( q = 0; q < n; ++q ) {
      cComplexType x = UL(p,q);
      UL(p,q) = UL(j,q);
      UL(j,q) = x;
    }
  }
}

