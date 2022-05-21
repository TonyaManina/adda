/* Inline complex functions, functions on length-3 real and complex vectors, and several auxiliary functions
 *
 * Copyright (C) ADDA contributors
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */
#ifndef __cmplx_h
#define __cmplx_h
#define _USE_MATH_DEFINES
// project headers
#include "const.h"    // for math constants
#include "types.h"    // for doublecomplex
// system headers
#include <math.h>
#include <string.h>   // for memcpy
#include <stdio.h>

// Uncomment this to turn off calculation of imExp using tables
//#define NO_IMEXP_TABLE

#ifdef USE_SSE3
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#endif

// useful macro for printing complex numbers and matrices
#define REIM(a) creal(a),cimag(a)
#define REIM3V(a) REIM((a)[0]),REIM((a)[1]),REIM((a)[2])

#ifndef NO_IMEXP_TABLE
void imExpTableInit(void);
doublecomplex imExpTable(double arg);
#endif
void imExp_arr(doublecomplex arg,int size,doublecomplex *c);

/* We do not use 'restrict' in the following functions since they are all inline - compiler will optimize the code
 * inside the calling function and decide whether the arrays can alias or not.
 */

extern doublecomplex Eye3[3][3];

//======================================================================================================================
// operations on complex numbers

static inline double cAbs2(const doublecomplex a)
// square of absolute value of complex number; |a|^2
{
	return creal(a)*creal(a) + cimag(a)*cimag(a);
}

//======================================================================================================================

static inline doublecomplex cSquare(const doublecomplex a)
{
	return creal(a)*creal(a)-cimag(a)*cimag(a)+2*creal(a)*cimag(a);
}

static inline doublecomplex cSqrtCut(const doublecomplex a)
// square root of complex number, with explicit handling of branch cut (not to depend on sign of zero of imaginary part)
/* It is designed for calculating normal component of the transmitted wavevector when passing through the plane
 * interface. However, such choice of branch cut (while physically correct) leads to all kind of weird consequences.
 *
 * For instance, the electric field above the interface for plane wave propagating from a slightly absorbing substrate
 * at large incident angle (larger than critical angle for purely real refractive index) is unexpectedly large. This
 * happens because the wave in the vacuum is inhomogeneous and the real part of wavevector is almost parallel to the
 * surface. So the field above the surface actually comes from distant points on the surface, which has much larger
 * amplitude of the incident wave from below (compared to that under the observation point). Since the distance along
 * the surface (or the corresponding slope) is inversely proportional to the imaginary part of the substrate refractive
 * index, the effect remains finite even in the limit of absorption going to zero. Therefore, in this case there exist
 * a discontinuity when switching from non-absorbing to absorbing substrate. Physically, this fact is a consequence of
 * the infinite lateral extent of the plane wave.
 *
 * Exactly the same issue exist when scattering into the absorbing medium is calculated. At large scattering angles the
 * amplitude becomes very large, which also amplifies a lot the calculated Csca.
 */
{
	if (cimag(a)==0) {
		if (creal(a)>=0) return sqrt(creal(a));
		else return I*sqrt(-creal(a));
	}
	else return csqrt(a);
}

//======================================================================================================================

static inline doublecomplex imExp(const double arg)
/* exponent of imaginary argument Exp(i*arg)
 * !!! should not be used in parameter parsing (table is initialized in VariablesInterconnect())
 */
{
#ifdef NO_IMEXP_TABLE
	/* We tried different standard options. (cos + I*sin) is almost twice slower than cexp, while sincos (GNU extension)
	 * is slightly faster (3.52 - 2.39 - 2.29 for matvec in test sparse runs, where about 1.23 is for non-exp part -
	 * median values over 10 runs). So we prefer to use standard cexp.
	 * When using table (below) the corresponding timing is 1.70.
	 */
	return cexp(I*arg);
#else
	return imExpTable(arg);
#endif
}

//======================================================================================================================
// operations on complex vectors

static inline void cvInit(doublecomplex a[static 3])
// set complex vector[3] to zero; a=0
{
	a[0] = 0;
	a[1] = 0;
	a[2] = 0;
}

//======================================================================================================================

static inline void vConj(const doublecomplex a[static 3],doublecomplex b[static 3])
// complex conjugate of vector[3]; b=a*
{
	b[0] = conj(a[0]);
	b[1] = conj(a[1]);
	b[2] = conj(a[2]);
}

//======================================================================================================================

static inline double cInvIm(const doublecomplex a)
// returns Im of inverse of a; designed to avoid under and overflows
{
	double tmp, a_RE, a_IM;
	a_RE=creal(a);
	a_IM=cimag(a);

	if (fabs(a_RE)>=fabs(a_IM)) {
		tmp=a_IM/a_RE;
		return (-tmp/(a_RE+a_IM*tmp));
	}
	else {
		tmp=a_RE/a_IM;
		return (-1/(a_RE*tmp+a_IM));
	}
}

//============================================================
static inline void vReal(const doublecomplex a[static 3],double b[static 3])
// takes real part of the complex vector; b=Re(a)
{
	b[0]=creal(a[0]);
	b[1]=creal(a[1]);
	b[2]=creal(a[2]);
}

//======================================================================================================================

static inline void cvBuildRe(const double a[static 3],doublecomplex b[static 3])
// builds complex vector from real part; b=a + i*0
{
	b[0]=a[0];
	b[1]=a[1];
	b[2]=a[2];
}

//======================================================================================================================

static inline void vInvRefl_cr(const double a[static 3],doublecomplex b[static 3])
/* reflects real vector with respect to the xy-plane and then inverts it, equivalent to reflection about the z-axis
 * result is stored into the complex vector
 */
{
	b[0]=-a[0];
	b[1]=-a[1];
	b[2]=a[2];
}

//======================================================================================================================

static inline void crCrossProd(const double a[static 3],const doublecomplex b[static 3],doublecomplex c[static 3])
// cross product of real and complex vector; c = a x b; !!! vectors must not alias
{
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
}

//======================================================================================================================

static inline void cvMultScal(const double a,const doublecomplex b[static 3],doublecomplex c[static 3])
// multiplication of vector[3] by real scalar; c=ab;
{
	c[0] = a*b[0];
	c[1] = a*b[1];
	c[2] = a*b[2];
}

//======================================================================================================================

static inline void cvMultScal_RVec(const doublecomplex a,const double b[static 3],doublecomplex c[static 3])
// complex scalar - real vector[3] multiplication; c=ab
{
	c[0] = a*b[0];
	c[1] = a*b[1];
	c[2] = a*b[2];
}

//======================================================================================================================

static inline void cvMultScal_cmplx(const doublecomplex a,const doublecomplex b[static 3],doublecomplex c[static 3])
// multiplication of vector[3] by complex scalar; c=ab
{
	c[0] = a*b[0];
	c[1] = a*b[1];
	c[2] = a*b[2];
}

//======================================================================================================================

static inline double cvNorm2(const doublecomplex a[static 3])
// square of the norm of a complex vector[3]
{
	return cAbs2(a[0]) + cAbs2(a[1]) + cAbs2(a[2]);
}


//======================================================================================================================

static inline doublecomplex cDotProd(const doublecomplex a[static 3],const doublecomplex b[static 3])
// conjugate dot product of two complex vector[3]; c=a.b = a[0]*b*[0]+...+a[2]*b*[2]; for one vector use cvNorm2
{
	return a[0]*conj(b[0]) + a[1]*conj(b[1]) + a[2]*conj(b[2]);
}

//======================================================================================================================

static inline double cDotProd_Im(const doublecomplex a[static 3],const doublecomplex b[static 3])
/* imaginary part of dot product of two complex vector[3]; c=Im(a.b)
 * It is not clear whether this way is faster than cimag(cDotProd)
 */
{
	return ( cimag(a[0])*creal(b[0]) - creal(a[0])*cimag(b[0])
	       + cimag(a[1])*creal(b[1]) - creal(a[1])*cimag(b[1])
	       + cimag(a[2])*creal(b[2]) - creal(a[2])*cimag(b[2]) );
}

//======================================================================================================================

static inline doublecomplex cDotProd_conj(const doublecomplex a[static 3],const doublecomplex b[static 3])
// dot product of two complex vector[3] without conjugation; a.(b*) = a[0]*b[0]+...+a[2]*b[2]
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

//======================================================================================================================

static inline void cvAdd(const doublecomplex a[static 3],const doublecomplex b[static 3],doublecomplex c[static 3])
// add two complex vector[3]; c=a+b;
{
	c[0] = a[0] + b[0];
	c[1] = a[1] + b[1];
	c[2] = a[2] + b[2];
}

//======================================================================================================================

static inline void cvSubtr(const doublecomplex a[static 3],const doublecomplex b[static 3],doublecomplex c[static 3])
// add two complex vector[3]; c=a-b;
{
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
}

//======================================================================================================================

static inline void cvAdd2Self(doublecomplex a[static 3],const doublecomplex b[static 3],const doublecomplex c[static 3])
// increment one complex vector[3] by sum of other two; a+=b+c
{
	a[0] += b[0] + c[0];
	a[1] += b[1] + c[1];
	a[2] += b[2] + c[2];
}

//======================================================================================================================

static inline doublecomplex crDotProd(const doublecomplex a[static 3],const double b[static 3])
// dot product of complex and real vectors[3]; a.b
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

//======================================================================================================================

static inline double crDotProd_Re(const doublecomplex a[static 3],const double b[static 3])
// real part of dot product of complex and real vectors[3]; Re(a.b)
{
	return creal(a[0])*b[0] + creal(a[1])*b[1] + creal(a[2])*b[2];
}

//======================================================================================================================

static inline void cvLinComb1(const doublecomplex a[static 3],const doublecomplex b[static 3],const double c1,
	doublecomplex c[static 3])
// linear combination of complex vectors[3] with real coefficient; second coefficient is unity; c=c1*a+b
{
	c[0] = c1*a[0] + b[0];
	c[1] = c1*a[1] + b[1];
	c[2] = c1*a[2] + b[2];
}

//============================================================

static inline void cvLinComb1_cmplx(doublecomplex a[static 3],doublecomplex b[static 3],const doublecomplex c1,
	doublecomplex c[static 3])
// linear combination of complex vectors[3] with complex coefficients; second coefficient is unity; c=c1*a+b
{
	c[0] = c1*a[0] + b[0];
	c[1] = c1*a[1] + b[1];
	c[2] = c1*a[2] + b[2];
}

//======================================================================================================================

static inline void cSymMatrVec(const doublecomplex matr[static 6],const doublecomplex vec[static 3],
	doublecomplex res[static 3])
// multiplication of complex symmetric matrix[6] by complex vec[3]; res=matr.vec; !!! vec and res must not alias
{
	res[0] = matr[0]*vec[0] + matr[1]*vec[1] + matr[2]*vec[2];
	res[1] = matr[1]*vec[0] + matr[3]*vec[1] + matr[4]*vec[2];
	res[2] = matr[2]*vec[0] + matr[4]*vec[1] + matr[5]*vec[2];
}

//======================================================================================================================

static inline void cSymMatrVecReal(const doublecomplex matr[static 6],const double vec[static 3],
	doublecomplex res[static 3])
// multiplication of complex symmetric matrix[6] by  vec[3]; res=matr.vec; !!! vec and res must not alias
{
	res[0] = matr[0]*vec[0] + matr[1]*vec[1] + matr[2]*vec[2];
	res[1] = matr[1]*vec[0] + matr[3]*vec[1] + matr[4]*vec[2];
	res[2] = matr[2]*vec[0] + matr[4]*vec[1] + matr[5]*vec[2];
}

//======================================================================================================================

static inline void cReflMatrVec(const doublecomplex matr[static 6],const doublecomplex vec[static 3],
	doublecomplex res[static 3])
/* multiplication of matrix[6] by complex vec[3]; res=matr.vec; passed components are the same as for symmetric matrix:
 * 11,12,13,22,23,33, but the matrix has the following symmetry - M21=M12, M31=-M13, M32=-M23
 * !!! vec and res must not alias
 */
{
	res[0] = matr[0]*vec[0] + matr[1]*vec[1] + matr[2]*vec[2];
	res[1] = matr[1]*vec[0] + matr[3]*vec[1] + matr[4]*vec[2];
	res[2] = - matr[2]*vec[0] - matr[4]*vec[1] + matr[5]*vec[2];
}

//======================================================================================================================

static inline void cReflMatrVecReal(const doublecomplex matr[static 6],const double vec[static 3],
	doublecomplex res[static 3])
/* multiplication of matrix[6] by complex vec[3]; res=matr.vec; passed components are the same as for symmetric matrix:
 * 11,12,13,22,23,33, but the matrix has the following symmetry - M21=M12, M31=-M13, M32=-M23
 * !!! vec and res must not alias
 */
{
	res[0] = matr[0]*vec[0] + matr[1]*vec[1] + matr[2]*vec[2];
	res[1] = matr[1]*vec[0] + matr[3]*vec[1] + matr[4]*vec[2];
	res[2] = - matr[2]*vec[0] - matr[4]*vec[1] + matr[5]*vec[2];
}

//======================================================================================================================
// operations on real vectors

static inline void vInit(double a[static 3])
// set real vector[3] to zero; a=0
{
	a[0] = 0;
	a[1] = 0;
	a[2] = 0;
}

//======================================================================================================================

static inline void vCopy(const double a[static 3],double b[static 3])
// copies one vector into another; b=a
{
	// can be rewritten through memcpy, but compiler should be able to do it itself if needed
	b[0]=a[0];
	b[1]=a[1];
	b[2]=a[2];
}

//======================================================================================================================

static inline void vAdd(const double a[static 3],const double b[static 3],double c[static 3])
// adds two real vectors; c=a+b
{
	c[0]=a[0]+b[0];
	c[1]=a[1]+b[1];
	c[2]=a[2]+b[2];
}

//======================================================================================================================

static inline void vSubtr(const double a[static 3],const double b[static 3],double c[static 3])
// subtracts two real vectors; c=a-b
{
	c[0]=a[0]-b[0];
	c[1]=a[1]-b[1];
	c[2]=a[2]-b[2];
}

//======================================================================================================================

static inline void vInvSign(double a[static 3])
// inverts the sign in the double vector[3]
{
	a[0]=-a[0];
	a[1]=-a[1];
	a[2]=-a[2];
}

//======================================================================================================================

static inline void vRefl(const double inc[static 3],double ref[static 3])
// reflects the incident vector 'inc' with respect to the xy-plane (inverts z-component)
{
	ref[0]=inc[0];
	ref[1]=inc[1];
	ref[2]=-inc[2];
}

//======================================================================================================================

static inline void vMultScal(const double a,const double b[static 3],double c[static 3])
// multiplication of real vector by scalar; c=a*b;
{
	c[0]=a*b[0];
	c[1]=a*b[1];
	c[2]=a*b[2];
}

//======================================================================================================================


static inline void vMultScalSelf(const double a,double b[static 3])
// multiplication of real vector by scalar; b*=a;
{
	b[0]*=a;
	b[1]*=a;
	b[2]*=a;
}

//============================================================
static inline void vMult(const double a[static 3],const double b[static 3],double c[static 3])
// multiplication of two vectors (by elements); c[i]=a[i]*b[i];
{
	c[0]=a[0]*b[0];
	c[1]=a[1]*b[1];
	c[2]=a[2]*b[2];
}

//======================================================================================================================

static inline bool vAlongZ(const double a[static 3])
// a robust (with respect to round-off errors) way to test that vector is along the z-axis (+ or -)
{
	return fabs(a[0])<ROUND_ERR && fabs(a[1])<ROUND_ERR;
}

//======================================================================================================================

static inline double DotProd(const double a[static 3],const double b[static 3])
// dot product of two real vectors[3]; use DotProd(x,x) to get squared norm
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

//======================================================================================================================

static inline double DotProdSquare(const double a[static 3],const double b[static 3])
// dot product of element-wise squares of two real vectors[3]
{
	return a[0]*a[0]*b[0]*b[0] + a[1]*a[1]*b[1]*b[1] + a[2]*a[2]*b[2]*b[2];
}

//======================================================================================================================

static inline double vNorm(const double a[static 3])
// norm of a real vector[3]
{
	return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

//============================================================

static inline double AngleCos(const double v[static 3], const double u[static 3])
// cosine of angle between vectors u and v
{
	return DotProd(v, u) / vNorm(v) * vNorm(u);
}

//======================================================================================================================

static inline void CrossProd(const double a[static 3],const double b[static 3],double c[static 3])
// cross product of two real vectors; c = a x b; !!! vectors must not alias
{
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
}

//======================================================================================================================
//2.01.22. static inline? static restrict?

static inline double AbsOutProd(const double a[static 3],const double b[static 3])
/* norm of outer product of two real vectors[3] (area of the corresponding parallelogram)
 * !!! pointers a and b must not alias !!! (does not make sense to be the same)
 */
{
	double x,y,z;
	x=a[1]*b[2]-a[2]*b[1];
	y=a[2]*b[0]-a[0]*b[2];
	z=a[0]*b[1]-a[1]*b[0];

	return sqrt(x*x+y*y+z*z);
}

//============================================================


static inline void vNormalize(double a[static 3])
// normalize real vector to have unit norm
{
	double c;
	c=1/sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
	a[0]*=c;
	a[1]*=c;
	a[2]*=c;
}

//======================================================================================================================

static inline void LinComb(const double a[static 3],const double b[static 3],const double c1,const double c2,
	double c[static 3])
// linear combination of real vectors[3]; c=c1*a+c2*b
{
	c[0]=c1*a[0]+c2*b[0];
	c[1]=c1*a[1]+c2*b[1];
	c[2]=c1*a[2]+c2*b[2];
}

//======================================================================================================================

static inline void OuterSym(const double a[static 3],double matr[static 6])
// outer product of real vector a with itself, stored in symmetric matrix matr
{
	matr[0] = a[0]*a[0];
	matr[1] = a[0]*a[1];
	matr[2] = a[0]*a[2];
	matr[3] = a[1]*a[1];
	matr[4] = a[1]*a[2];
	matr[5] = a[2]*a[2];
}

//======================================================================================================================

static inline double TrSym(const double a[static 6])
// trace of a symmetric matrix stored as a vector of size 6
{
	return (a[0]+a[2]+a[5]);
}

//======================================================================================================================

static inline double QuadForm(const double matr[static 6],const double vec[static 3])
// value of a quadratic form matr (symmetric matrix stored as a vector of size 6) over a vector vec;
{
	return ( vec[0]*vec[0]*matr[0] + vec[1]*vec[1]*matr[2] + vec[2]*vec[2]*matr[5]
	       + 2*(vec[0]*vec[1]*matr[1] + vec[0]*vec[2]*matr[3] + vec[1]*vec[2]*matr[4]) );
}

//======================================================================================================================

static inline void MatrVec(double matr[static 3][3],const double vec[static 3],double res[static 3])
// multiplication of matrix[3][3] by vec[3] (all real); res=matr.vec;
{
	res[0]=matr[0][0]*vec[0]+matr[0][1]*vec[1]+matr[0][2]*vec[2];
	res[1]=matr[1][0]*vec[0]+matr[1][1]*vec[1]+matr[1][2]*vec[2];
	res[2]=matr[2][0]*vec[0]+matr[2][1]*vec[1]+matr[2][2]*vec[2];
}

//======================================================================================================================

static inline void MatrVecComplex(doublecomplex matr[static 3][3],const doublecomplex vec[static 3],doublecomplex res[static 3])
// multiplication of matrix[3][3] by vec[3] (all real); res=matr.vec;
{
	res[0]=matr[0][0]*vec[0]+matr[1][0]*vec[1]+matr[2][0]*vec[2];
	res[1]=matr[0][1]*vec[0]+matr[1][1]*vec[1]+matr[2][1]*vec[2];
	res[2]=matr[0][2]*vec[0]+matr[1][2]*vec[1]+matr[2][2]*vec[2];
}

//======================================================================================================================

static inline void MatrColumn(double matr[static 3][3],const int ind,double vec[static 3])
// get ind's column of matrix[3][3] and store it into vec[3] (all real, ind starts from zero); vec=matr[.][ind];
{
	vec[0]=matr[0][ind];
	vec[1]=matr[1][ind];
	vec[2]=matr[2][ind];
}

//======================================================================================================================

static inline void MatrPlainTo3x3(doublecomplex plain[static 9], doublecomplex m[static 3][3]) {
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			m[i][j] = plain[j*3 + i];
}

//======================================================================================================================

static inline void DblPlainToCmplx3x3(double plain[static 9], doublecomplex m[static 3][3]) {
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			m[i][j] = plain[j*3 + i];
}


//======================================================================================================================

static inline void MatrDotProd(int n, doublecomplex A[][n], doublecomplex B[][n], doublecomplex C[][n]) {
	//C=A*B
	for (int j = 0; j < n; j++) // A * x_j = e_j
		for (int i = 0; i < n; i++) {
			C[i][j]=0;
			for (int k = 0; k < n; k++){
				C[i][j] += A[k][j] * B[i][k];
			}
		}
}

//======================================================================================================================

static inline void MatrSet(doublecomplex m[static 3][3], doublecomplex val) {
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			m[i][j] = val;
}

//======================================================================================================================

static inline void MatrVecMul(int n, doublecomplex A[][n], doublecomplex vec[3],doublecomplex res[3]) {
	for (int row = 0; row < n; row++) {
		res[row] = 0;
		for (int col = 0; col < n; col++){
			res[row] += A[col][row]*vec[col];
		}
	}
}

//======================================================================================================================

static inline void DblMatrDotProd(int n, double A[][n], double B[][n], double C[][n]) {
	//C=A*B
	double a,b,c;
	for (int row = 0; row < n; row++) { // A * x_j = e_j
		for (int col = 0; col < n; col++) {
			C[col][row]=0;
			for (int k = 0; k < n; k++){
				a = A[k][row];
				b = B[col][k];
				c = a * b;
				C[col][row] += c;
			}
		}
	}
}


//======================================================================================================================

static inline void DblMatrSet(double m[static 3][3], double val) {
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			m[i][j] = val;
}

//======================================================================================================================
static inline void MatrMul(doublecomplex m[static 3][3], doublecomplex val) {
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			m[i][j] *= val;
}

//======================================================================================================================

static inline doublecomplex MatrDet(doublecomplex m[static 3][3]) {
	return
			m[0][0]*m[1][1]*m[2][2] +
			m[1][0]*m[2][1]*m[0][2] +
			m[2][0]*m[0][1]*m[1][2] -
			m[2][0]*m[1][1]*m[0][2] -
			m[1][0]*m[0][1]*m[2][2] -
			m[0][0]*m[1][2]*m[2][1];
}

//======================================================================================================================

static inline void CubicSolver(doublecomplex a, doublecomplex b, doublecomplex c, doublecomplex d, doublecomplex root[static 3]) {
	doublecomplex p = c/a - b*b/(3*a*a);
	doublecomplex q = 2*b*b*b/(27*a*a*a)-b*c/(3*a*a)+d/a;
	doublecomplex tmp = csqrt(q*q/4 + p*p*p/27);
	doublecomplex exponenta = cexp(2*I*PI/3);
	doublecomplex z1 = cpow(-q/2 - tmp, 1.0 / 3.0);
	doublecomplex z2 = z1*exponenta, z3 = z2*exponenta;
	doublecomplex y1 = z1 - p/(3*z1), y2 = z2 - p/(3*z2), y3 = z3 - p/(3*z3);
	root[0] = y1-b/(3*a); root[1] = y2-b/(3*a); root[2] = y3-b/(3*a);
}

//======================================================================================================================

static inline void Eigenvalues(doublecomplex m[static 3][3], doublecomplex lambda[static 3]) {
	doublecomplex Det = MatrDet(m);
	doublecomplex Tr = m[0][0]+m[1][1]+m[2][2];
	doublecomplex c = m[0][0]*m[1][1]-m[0][1]*m[1][0]+m[1][1]*m[2][2]-m[1][2]*m[2][1]+m[0][0]*m[2][2]-m[0][2]*m[2][0];
	CubicSolver(-1.0, Tr, -c, Det, lambda);
}

//======================================================================================================================

static inline void MatrCopy(int n, doublecomplex dest[][n], doublecomplex src[][n]) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			dest[i][j] = src[i][j];
}

//======================================================================================================================

static inline void DblMatrCopy(int n, double dest[][n], double src[][n]) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			dest[i][j] = src[i][j];
}

//======================================================================================================================

static inline void MatrInverse(doublecomplex matr[static 3][3], doublecomplex inv[static 3][3])
// get ind's column of matrix[3][3] and store it into vec[3] (all real, ind starts from zero); vec=matr[.][ind];
{
	doublecomplex det = MatrDet(matr);
	doublecomplex deltaMatr[3][3];

	for (int j = 0; j < 3; j++) // A * x_j = e_j
		for (int i = 0; i < 3; i++) { // x_j [i]
			MatrCopy(3, deltaMatr, matr);
			for (int k = 0; k < 3; k++)
				deltaMatr[i][k] = (k == j) ? 1.0 : 0.0;

			doublecomplex delta = MatrDet(deltaMatr);
			for (int k = 0; k < 3; k++)
				inv[j][i] = delta / det;

		}
}

//======================================================================================================================

static inline void MatrAdd(doublecomplex dest[static 3][3], doublecomplex src[static 3][3]) {
	for (int j = 0; j < 3; j++) // A * x_j = e_j
		for (int i = 0; i < 3; i++) // x_j [i]
			dest[i][j] += src[i][j];
}

//======================================================================================================================
static inline void MatrSubtract(doublecomplex dest[static 3][3], doublecomplex src[static 3][3]) {
	for (int j = 0; j < 3; j++) // A * x_j = e_j
		for (int i = 0; i < 3; i++) // x_j [i]
			dest[i][j] -= src[i][j];
}
//======================================================================================================================
static inline void MatrTranspose(doublecomplex dest[static 3][3], doublecomplex src[static 3][3]) {
	for (int j = 0; j < 3; j++) // A * x_j = e_j
		for (int i = 0; i < 3; i++) // x_j [i]
			dest[i][j] = src[j][i];
}
//======================================================================================================================
static inline void DblMatrTranspose(int n, double dest[][n], double src[][n]) {
	for (int j = 0; j < n; j++) // A * x_j = e_j
		for (int i = 0; i < n; i++) // x_j [i]
			dest[i][j] = src[j][i];
}
//======================================================================================================================
static inline void DblPrintMatr(int n, double v[][n]) {
	double d = 0;
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			d = v[i][j];
			printf("%.3f ", d);
		}
		printf("\n");
	}
}
//======================================================================================================================
static inline void DblcomplexPrintMatr(int n, doublecomplex v[][n]) {
	doublecomplex d = 0;
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			d = v[i][j];
			printf("%.3f + i%.3f ", creal(d), cimag(d));
		}
		printf("\n");
	}
}
//======================================================================================================================
static inline void MatrixRoot(doublecomplex M[static 3][3], doublecomplex lambda[static 3], doublecomplex Mroot[static 3][3]) {
	doublecomplex E1[3][3] = {
			{lambda[0],0,0},
			{0,lambda[0],0},
			{0,0,lambda[0]}
	};
	doublecomplex E2[3][3] = {
			{lambda[1],0,0},
			{0,lambda[1],0},
			{0,0,lambda[1]}
	};
	doublecomplex E3[3][3] = {
			{lambda[2],0,0},
			{0,lambda[2],0},
			{0,0,lambda[2]}
	};
	doublecomplex boof1[3][3], boof2[3][3], boof3[3][3];
	MatrCopy(3, boof1, M); MatrSubtract(boof1, E1);
	MatrCopy(3, boof2, M); MatrSubtract(boof2, E2);
	MatrCopy(3, boof3, M); MatrSubtract(boof3, E3);
	doublecomplex C[3][3];
	MatrDotProd(3, boof2,boof3,C);
	doublecomplex c = 1.0/((lambda[0]-lambda[1])*(lambda[0]-lambda[2]));
	MatrMul(C,c*csqrt(lambda[0]));
	MatrAdd(Mroot,C);
	MatrDotProd(3, boof1,boof3,C);
	c = 1.0/((lambda[1]-lambda[2])*(lambda[1]-lambda[0]));
	MatrMul(C,c*csqrt(lambda[1]));
	MatrAdd(Mroot, C);
	MatrDotProd(3, boof1,boof2,C);
	c = 1.0/((lambda[2]-lambda[0])*(lambda[2]-lambda[1]));
	MatrMul(C,c*csqrt(lambda[2]));
	MatrAdd(Mroot, C);
}

//======================================================================================================================
static inline void Permutate(double vec[static 3],const int ord[static 3])
// permutate double vector vec using permutation ord
{
	double buf[3];

	memcpy(buf,vec,3*sizeof(double));
	vec[0]=buf[ord[0]];
	vec[1]=buf[ord[1]];
	vec[2]=buf[ord[2]];
}

//======================================================================================================================

static inline void Permutate_i(int vec[static 3],const int ord[static 3])
// permutate int vector vec using permutation ord
{
	int buf[3];

	memcpy(buf,vec,3*sizeof(int));
	vec[0]=buf[ord[0]];
	vec[1]=buf[ord[1]];
	vec[2]=buf[ord[2]];
}

//======================================================================================================================
// Auxiliary functions

static inline double Deg2Rad(const double deg)
// transforms angle in degrees to radians
{
	return (PI_OVER_180*deg);
}

//======================================================================================================================

static inline double Rad2Deg(const double rad)
// transforms angle in radians to degrees
{
	return (INV_PI_180*rad);
}

//======================================================================================================================

static inline bool TestBelowDeg(const double deg)
/* tests if the direction is below the substrate using the degree theta in degrees;
 * if unsure (within rounded error) returns false (above)
 */
{
	return fabs(fmod(fabs(deg),360)-180) < 90*(1-ROUND_ERR);
}

//======================================================================================================================
// functions used for substrate

//======================================================================================================================

static inline doublecomplex FresnelRS(const doublecomplex ki,const doublecomplex kt)
/* reflection coefficient for s-polarized wave (E perpendicular to the main plane),
 * ki,kt are normal (positive) components of incident and transmitted wavevector (with arbitrary mutual scaling)
 */
{
	return (ki-kt)/(ki+kt);
}

//======================================================================================================================

static inline doublecomplex FresnelTS(const doublecomplex ki,const doublecomplex kt)
/* transmission coefficient for s-polarized wave (E perpendicular to the main plane),
 * ki,kt are normal (positive) components of incident and transmitted wavevector (with arbitrary mutual scaling)
 */
{
	return 2*ki/(ki+kt);
}

//======================================================================================================================

static inline doublecomplex FresnelRP(const doublecomplex ki,const doublecomplex kt,const doublecomplex mr)
/* reflection coefficient for p-polarized wave (E parallel to the main plane),
 * ki,kt are normal (positive) components of incident and transmitted wavevector (with arbitrary mutual scaling)
 * mr is the ratio of refractive indices (mt/mi)
 */
{
	return (mr*mr*ki-kt)/(mr*mr*ki+kt);
}

//======================================================================================================================

static inline doublecomplex FresnelTP(const doublecomplex ki,const doublecomplex kt,const doublecomplex mr)
/* transmission coefficient for p-polarized wave (E parallel to the main plane),
 * ki,kt are normal (positive) components of incident and transmitted wavevector (with arbitrary mutual scaling)
 * mr is the ratio of refractive indices (mt/mi)
 */
{
	return 2*mr*ki/(mr*mr*ki+kt);
}

#ifdef USE_SSE3

//======================================================================================================================

static inline __m128d cmul(__m128d a,__m128d b)
// complex multiplication
{
	__m128d t;
	t = _mm_movedup_pd(a);
	a = _mm_shuffle_pd(a,a,3);
	t = _mm_mul_pd(b,t);
	b = _mm_shuffle_pd(b,b,1);
	b = _mm_mul_pd(a,b);
	return _mm_addsub_pd(t,b);
}

//======================================================================================================================

static inline __m128d cadd(__m128d a,__m128d b)
// complex addition
{
	return _mm_add_pd(a,b);
}

#endif // USE_SSE3

#endif // __cmplx_h
