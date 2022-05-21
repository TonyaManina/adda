/* All the initialization is done here before actually calculating internal fields,
 * includes calculation of couple constants
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
#define _USE_MATH_DEFINES
#include "const.h" // keep this first
// project headers
#include "cmplx.h"
#include "comm.h"
#include "crosssec.h"
#include "debug.h"
#include "fft.h"
#include "interaction.h"
#include "io.h"
#include "memory.h"
#include "oclcore.h"
#include "Romberg.h"
#include "timing.h"
#include "vars.h"
#include "volfrac.h"
#include "qr.h"
#include "tagaki_factor.h"

// system headers
#include <math.h>
#include <stdlib.h>
#include <string.h>

//======================================================================================================================

const double epsEq = 10e-12;
const double cubeOrderedPoints[8][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,0,1},{1,1,0},{0,1,1},{1,1,1}};
const int cubeOrderedEdges[6][4] = {{0,2,5,1}, {0,3,6,2}, {3,4,7,6}, {4,1,5,7}, {0,1,4,3}, {2,6,7,5}};
const double CubeCenter[3] = {0.5, 0.5, 0.5};
const double CubeEdgeNorm[6][3] = {{0,0,-1},{-1,0,0},{0,0,1},{1,0,0},{0,-1,0},{0,1,0}};
const double check = 0.00000001;

//======================================================================================================================

void PrintVector(const double p[static 3]) {
	printf("(%.3f, %.3f, %.3f)", p[0], p[1], p[2]);
}

//======================================================================================================================

bool vEq(double v[3], double u[3]) {
	double diff[3];
	vSubtr(v, u, diff);
	return DotProd(diff, diff) < epsEq;
}

//==========================================================

void BubbleSort(int* indxs, double* values, int n) {
	double boof;
	for (int i = 0; i<n; i++){
		for (int j=i+1; j<n; j++){
			double vi = values[indxs[i]], vj = values[indxs[j]];
			if (vi<vj) {
				boof=indxs[i];
				indxs[i]=indxs[j];
				indxs[j]=boof;
			}
		}
	}
}

//======================================================================================================================

bool PointUpperPlane(const double p[static 3], const double plane_n[static 3]) {
	return DotProd(p, plane_n) > 1.0;
}

//======================================================================================================================

void FillCubePlaneIntHelperArray(double nv[8], double n[3]) {
	for (int i = 0; i < 8; i++) {
		nv[i] = DotProd(cubeOrderedPoints[i], n);
	}
}

//======================================================================================================================

int CubePlaneRawIntersections(const double n[3], double intersections[24][3], int edgeCounts[6]) {
	int i_idx = 0;
	double nv[8];

	FillCubePlaneIntHelperArray(nv, n);
	for (int nEdge = 0; nEdge < 6; nEdge++) {
		int edgeIntCount = 0;
		for (int j = 0; j < 4; j++) {
			int k = (j+1) % 4;
			int v1 = cubeOrderedEdges[nEdge][j];
			int v2 = cubeOrderedEdges[nEdge][k];
			if (EdgeIn(v1,v2,nv,intersections[i_idx])) {
				i_idx++;
				edgeIntCount++;
			}

		}

		edgeCounts[nEdge] = edgeIntCount;
	}

	return i_idx;
}

//======================================================================================================================

void PointsCenter(const double p[24][3], int n, double c[3]) {
	vInit(c);
	for (int i = 0; i < n; i++)
		vAdd(c, p[i], c);

	vMultScal(1.0 / (double) n, c, c);
}

//======================================================================================================================

double PointWalkOrd(const double u0[3], const double u[3], const double n[3]) {
	double vecProd[3];

	if (vEq(u0, u))
		return 2;

	CrossProd(u0, u, vecProd);
	double s = DotProd(vecProd, n);
	if (s > 0) s = 1;
	if (s < 0) s = -1;

	return s * (1.0 + AngleCos(u0, u));
}

//======================================================================================================================

int ReorderPoints(const double p[24][3], int k, const double n[3], const double pReordered[24][3]) {
	double c[3] = {0,0,0};
	double u[24][3];
	double ords[24];
	int walkIndexes[24];
	if (k == 0) return 0;

	// switch to figure center coords
	PointsCenter(p, k, c);
	for (int i = 0; i < k; i++)
		vSubtr(p[i], c, u[i]);

	// calc ords for all vectors
	ords[0] = 2.0; walkIndexes[0] = 0;
	for (int i = 1; i < k; i++) {
		ords[i] = PointWalkOrd(u[0], u[i], n);
		walkIndexes[i] = i;
	}

	// reorder based on ords[], descending
	BubbleSort(walkIndexes, ords, k);

	// remove duplicate points and place to pReordered
	int prevApproved = 0, countPlaced = 1;
	vCopy(u[walkIndexes[0]], pReordered[0]);

	for (int i = 1; i < k; i++) {
		int curIndex = walkIndexes[i];
		if (vEq(u[curIndex], u[prevApproved]) != true) {
			vCopy(u[walkIndexes[i]], pReordered[countPlaced]);
			prevApproved = curIndex;
			countPlaced++;
		}
	}

	// switch back to 0,0,0 coords
	for (int i = 0; i < countPlaced; i++)
		vAdd(pReordered[i], c, pReordered[i]);

	return countPlaced;
}

//======================================================================================================================

void CalculationOfLsTensor(const double p[24][3], int k, const double n[3], double L[9] ) {
	if (k<3) return;
	//calculation of h
	double hlog[3] = {0,0,0}, q[3]= {0,0,0}, c[3], boof[3]={0,0,0};
	double logarifm;
	double u[24][3];

	c[0]=0.5; c[1]=0.5; c[2]=0.5;
    for (int j=0; j<k; j++)
    	 vSubtr(p[j],c,u[j]);
    vSubtr(p[0],c,u[k]);

    for (int j = 0; j<k; j++){
    	vSubtr(u[j+1], u[j], boof);
    	double boofNorm = vNorm(boof);
    	vMultScal(1.0 / boofNorm, boof, q);
    	logarifm=log((DotProd(u[j+1],q)+vNorm(u[j+1]))/(DotProd(u[j],q)+vNorm(u[j])));
    	vMultScal(logarifm,q,q);
    	vAdd(hlog,q,hlog);
    }

	double omega = 0; //calculation of omega
	for (int j=0; j<k; j++) {
		//vSubtr(p[j],c,u[j]);
    	double uNorm=vNorm(u[j]);
    	vMultScal(1.0 / uNorm, u[j], u[j]);
    }
    for (int j=0; j<k-2; j++) {
    	double prod[3];
    	CrossProd(u[j+1], u[j+2], prod);
    	double f = DotProd(u[j], prod);
    	double g = 1 + DotProd(u[j],u[j+1]) + DotProd(u[j],u[j+2]) + DotProd(u[j+1],u[j+2]);

    	if (fabs(g)<check && f>=0) omega += M_PI;
    	if (fabs(g)<check && f<0) omega -= M_PI;
    	if (fabs(g)>=check) {
    		if (f>=0) {
    			if (g>=0) {
    				omega += 2*atan(f/g);
    			}
    			else {
    				omega += 2*(atan(f/g)+M_PI);
    				}
    		}
    		else {
    			if (g>=0) {
    				omega += 2*atan(f/g);
    			}
    			else {
    				omega += 2*(atan(f/g)-M_PI);
    				}
    		}
    	}
    }

    double prod[3], omegan[3]; //Calculation of Ls
    CrossProd(n,hlog,prod);
    vMultScal(omega, n, omegan);
    vAdd(prod, omegan, prod);

    L[0]+=prod[0]*n[0]; L[1]+=prod[0]*n[1]; L[2]+=prod[0]*n[2];
    L[3]+=prod[1]*n[0]; L[4]+=prod[1]*n[1]; L[5]+=prod[1]*n[2];
    L[6]+=prod[2]*n[0]; L[7]+=prod[2]*n[1]; L[8]+=prod[2]*n[2];
}
//======================================================================================================================
void PolarizabilityCalc(doublecomplex refind, double vf, doublecomplex alpha[3][3], double n[3] ) {
	printf("refind = %.5f+I%.5f; vf = %.5f\n", creal(refind), cimag(refind), vf);
	double Ls[9] = {0,0,0,0,0,0,0,0,0};
	doublecomplex T[3][3], chi_eff[3][3], LsMatr[3][3]; //chi of outer space
	double UnsortedEdgePoints[6][24][3], SortedEdgePoints[6][24][3];
	double intersections[24][3], pOrdered[24][3], pOrdered_inv[24][3];
	int counting = 0, accumCounting = 0, NumberOfPoints[6], edgeIntCounts[6]; // separate counts of intersections (by edges)
	bool CenterUpperPlane = PointUpperPlane(CubeCenter, n);

	int k = CubePlaneRawIntersections(n, intersections, edgeIntCounts);
	k = ReorderPoints(intersections, k, n, pOrdered);
	for (int j=0; j<6; j++){
		counting=0;
		for (int i=0; i<4; i++){
			double* a;
			a = cubeOrderedPoints[cubeOrderedEdges[j][i]];
			if (PointUpperPlane(a,n) != CenterUpperPlane){
				vCopy(a, UnsortedEdgePoints[j][counting++]);
			}
		}

		for (int l=0; l<edgeIntCounts[j]; l++){
			vCopy(intersections[accumCounting + l], UnsortedEdgePoints[j][counting + l]);
		}

		counting += edgeIntCounts[j];
		accumCounting += edgeIntCounts[j];
		NumberOfPoints[j]=counting;
	}

	for (int i=0; i<6; i++)
		NumberOfPoints[i] = ReorderPoints(UnsortedEdgePoints[i], NumberOfPoints[i], CubeEdgeNorm[i], SortedEdgePoints[i]);

	for (int i=0; i<6; i++)
		CalculationOfLsTensor(SortedEdgePoints[i], NumberOfPoints[i], CubeEdgeNorm[i], Ls );

	if (CenterUpperPlane == false){
		vInvSign(n);
		for (int i=0; i<k; i++) vCopy(pOrdered[i],pOrdered_inv[k-i-1]);
		CalculationOfLsTensor(pOrdered_inv, k, n, Ls );
	} else {
		CalculationOfLsTensor(pOrdered, k, n, Ls );
	}
	vNormalize(n);
	//DblPlainToCmplx3x3(Ls,LsMatr);
	LsMatr[0][0] = Ls[0]; LsMatr[1][0] = Ls[1]; LsMatr[2][0] = Ls[2];
	LsMatr[0][1] = Ls[3]; LsMatr[1][1] = Ls[4]; LsMatr[2][1] = Ls[5];
	LsMatr[0][2] = Ls[6]; LsMatr[1][2] = Ls[7]; LsMatr[2][2] = Ls[8];
	doublecomplex LpMatr[3][3];
	MatrCopy(3,LpMatr,Eye3); MatrMul(LpMatr,4*M_PI/3.0);
	MatrSubtract(LpMatr,LsMatr);
	//MatrCopy(3,LsMatr,Eye3); MatrMul(LsMatr,4*M_PI/3.0); //сделала матрицу Ls для теста 4pi/3 * I
	doublecomplex temp[3][3], temp1[3][3], alphaT[3][3];
	doublecomplex M=(SO_B1*kd*kd+I*2*kd*kd*kd/3)*vf;
	//doublecomplex M = (0.278405 + I*0.048897524);
	doublecomplex chi_p = (refind*refind - 1.0)/ (4.0 * M_PI); // 1.0017 refind (как в стандартном запуске)
	//MatrSet(T, 0);
	//for (int i = 0; i<3; i++)
	//	for (int j=0; j<3; j++)
	//		T[j][i] = n[i]*n[j]; //заполняю Т, i - строка j - столбец
	//MatrMul(T, 1.0 / (refind * refind) - 1.0); //умножаю Т на 1/(e-1) (в числителе e среды = 1 поэтому 1/e)
//	MatrAdd(T, Eye3);// T = E + n n / (refind*refind-1)
	doublecomplex alphamaster;
	alphamaster = (3*dipvol*(refind*refind - 1)/(4*M_PI*(refind*refind+2)))/(1-3*(refind*refind-1)*M/(4*M_PI*(refind*refind+2)));
	MatrCopy(3, chi_eff, Eye3);
	MatrMul(chi_eff, vf * chi_p); //chi_eff[i*3+j]=(1-vf)*chi_out[i]*(i==j ? 1.0 : 0.0)+vf*1*T[i*3+j]; m_hoff not considered yet
	MatrCopy(3, temp, Eye3);
	MatrMul(temp,-M*vf);// Создала матрицу -М*I = temp
	MatrAdd(temp,LpMatr); // (Ls-MI)
	MatrMul(temp,chi_p); //(Ls-MI)*chi_s
	//MatrDotProd(3,temp,T,temp1); //temp1 = ((Ls-MI)*chi_s)* T
	MatrAdd(temp, Eye3); // temp1 = ((Ls-MI)*chi_s)* T + I (остальное не учитывается, потому что chi_p = 0
	MatrInverse(temp,temp1); //temp = (temp1)^-1
	doublecomplex test[3][3];
	MatrDotProd(3,temp1,temp,test);
	MatrDotProd(3, chi_eff,temp1,alpha); // alpha = chi_eff * temp
	MatrMul(alpha,dipvol); //alpha = alpha * dipvol
    MatrTranspose(alphaT,alpha);
    MatrAdd(alpha,alphaT);
    MatrMul(alpha,1.0/2.0);
}
void TestPolCalc() {
	double n[3] = {1,2,3};
	doublecomplex m = 1.0017, alpha[3][3];
	double vf = 1;
	PolarizabilityCalc(m, vf, alpha, n );
	int sdfsfs;

}
void Testmatrinv() {
	doublecomplex A[3][3] = {
			{3,98,1},
			{7,6,12},
			{4,51,8}
	};
	doublecomplex B[3][3], C[3][3];
	MatrInverse(A,B);
    MatrDotProd(3,A,B,C);
    int dfdfsd = 8;
}





