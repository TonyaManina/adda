/*
 * volfrac.h
 *
 *  Created on: 20 апр. 2022 г.
 *      Author: Workstation
 */

#ifndef SRC_VOLFRAC_H_
#define SRC_VOLFRAC_H_

bool PointUpperPlane(const double p[static 3], const double plane_n[static 3]);
int CountEdgePlaneIntersections(int edgeIdx, const double n[3], const double intersections[2][3]);
void TakagiDecomposition(doublecomplex alpha[3][3], doublecomplex beta[3][3], doublecomplex betaT[3][3]);
void PolarizabilityCalc(doublecomplex refind, double vf, doublecomplex alpha[3][3], double n[3]);
bool EdgeIn(int i,int j,double nv[8],double res[3]);
void TestPolCalc(void);
void Testmatrinv(void);

#endif /* SRC_VOLFRAC_H_ */
