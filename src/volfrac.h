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
void TestVolFrac(void);

bool EdgeIn(int i,int j,double nv[8],double res[3]);


#endif /* SRC_VOLFRAC_H_ */
