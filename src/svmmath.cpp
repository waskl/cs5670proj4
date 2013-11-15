/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #1:
 * svmmath.cpp
 *		a routine for intersecting >2 lines (for vanishing point
 *		computation);
 *		routines for computing the homography for the reference
 *		plane and arbitrary polygons
 **************************************************************/

#pragma warning(disable : 4996)

#include "Eigen/Core"
#include "MinEig.h"

#include "svmmath.h"
#include "jacob.h"
#include "vec.h"
#include <cstring>
#include <cstdio>
#include <assert.h>
#include <iostream>

using namespace Eigen;
using namespace std;

//
// TODO 1: BestFitIntersect()
//		Given lines, the list of 3 or more lines to be intersected,
//		find the best fit intersection point.
//		See http://www-2.cs.cmu.edu/~ph/869/www/notes/vanishing.txt.
//	
SVMPoint BestFitIntersect(const std::list<SVMLine> &lines, int imgWidth, int imgHeight)
{
    // check
    if (lines.size() < 2)
	{
            fprintf(stderr, "Not enough lines to compute the best fit.");
            abort();
	}

    SVMPoint bestfit;
    list<SVMLine>::const_iterator iter;

    // To accumulate stuff
    typedef Matrix<double, Dynamic, 3, RowMajor> Matrix3;

    int numLines = (int) lines.size();
	Matrix3 A = Matrix3::Zero(3, 3);	
    //Matrix3 A = Matrix3::Zero(numLines, 3);	
	
	//iter=lines.begin();
	//cout<<"Iter pt1,2 "<<iter->pnt1<<"\t"<<iter->pnt2<<"\n\n";
	
    // Transformation for numerical stability

    // Note: iterate through the lines list as follows:
    //		for (iter = lines.begin(); iter != lines.end(); iter++) {
    //			...iter is the pointer to the current line...
    //		}
    // Note: Function to find eigenvector with smallest eigenvalue is MinEig(A, eval, evec)
    //
    /******** BEGIN TODO ********/

	for (iter = lines.begin(); iter != lines.end(); iter++) {
    //			...iter is the pointer to the current line...

		SVMPoint p1, p2;
		p1.u = iter->pnt1->u;	p2.u = iter->pnt2->u; 
		p1.v = iter->pnt1->v;	p2.v = iter->pnt2->v;
		p1.w = iter->pnt1->w;	p2.w = iter->pnt2->w;

		p1.u-=imgWidth/2;				p2.u-=imgWidth/2;
		p1.v-=imgHeight/2;				p2.v-=imgHeight/2;
		p1.w=(imgWidth+imgHeight)/4;	p2.w=(imgWidth+imgHeight)/4;

		Vec3<double> Cross, P1, P2;
		P1[0] = p1.u;	P2[0] = p2.u;
		P1[1] = p1.v;	P2[1] = p2.v;
		P1[2] = p1.w;	P2[2] = p2.w;

		//cout<<"P1 "<<P1[0]<<" "<<P1[1]<<" "<<P1[2]<<"\n";
		//cout<<"P2 "<<P2[0]<<" "<<P2[1]<<" "<<P2[2]<<"\n";

		Cross = cross(P1,P2);

		//cout<<"Cross "<<Cross[0]<<" "<<Cross[1]<<" "<<Cross[2]<<"\n";

//		Cross[0] = p1.v*p2.w - p1.w*p2.v;
//		Cross[1] = p1.w*p2.u - p1.u*p2.w;
//		Cross[2] = p1.u*p2.v - p1.v*p2.u;
		
		A(0,0) += Cross[0]*Cross[0];	A(0,1) += Cross[0]*Cross[1];	A(0,2) += Cross[0]*Cross[2];
		A(1,0) += Cross[0]*Cross[1];	A(1,1) += Cross[1]*Cross[1];	A(1,2) += Cross[1]*Cross[2];
		A(2,0) += Cross[0]*Cross[2];	A(2,1) += Cross[1]*Cross[2];	A(2,2) += Cross[2]*Cross[2];
		

		//double Mag = sqrt( p1.u*p1.u + p1.v*p1.v + 1.0 );
		//p1.u/=Mag;	p1.v/=Mag;	p1.w/=Mag;

		//Mag = sqrt( p2.u*p2.u + p2.v*p2.v + 1.0 );
		//p2.u/=Mag;	p2.v/=Mag;	p2.w/=Mag;



		//cout<<"Iter pt1,2 "<<iter->pnt1->u<<","<<iter->pnt1->v<<","<<iter->pnt1->w<<"\t"<<iter->pnt2->u<<","<<iter->pnt2->v<<","<<iter->pnt2->w<<"\n\n";
    }

	double* evec = new double[3];
	double eval = 0;

	MinEig(A, eval, evec);

printf("TODO: %s:%d\n", __FILE__, __LINE__); 

    /******** END TODO ********/
	bestfit.u = evec[0] / evec[2];
	bestfit.v = evec[1] / evec[2];
	bestfit.w = 1.0;
	bestfit.u = bestfit.u*(imgWidth + imgHeight)/4 + imgWidth/2;
	bestfit.v = bestfit.v*(imgWidth + imgHeight)/4 + imgHeight/2;

    return bestfit;
}


//
// TODO 2: ConvertToPlaneCoordinate()
//		Given a plane defined by points, converts their coordinates into
//		plane coordinates. See the following document for more detail.
//		http://www.cs.cornell.edu/courses/cs4670/2012fa/projects/p4/homography.pdf.
//      The final divisors you apply to the u and v coordinates should be saved uScale and vScale
//
void ConvertToPlaneCoordinate(const vector<SVMPoint>& points, vector<Vec3d>& basisPts, double &uScale, double &vScale)
{
    int numPoints = points.size();
	
	for (int i=0; i < numPoints; i++) {
            Vec3d tmp = Vec3d(points[i].X, points[i].Y, points[i].Z);
            basisPts.push_back(tmp);
    }

	Vec3d Ex;
	Ex = basisPts[1]-basisPts[0];

	int Q = 2;
	double CloseToZero = 1;

	for (int i=2; i < numPoints; i++) {
		
        double Dot = (basisPts[i].operator[](0) - basisPts[0].operator[](0)) * Ex.operator[](0);
		Dot +=(basisPts[i].operator[](1) - basisPts[0].operator[](1)) * Ex.operator[](1);
		Dot +=(basisPts[i].operator[](2) - basisPts[0].operator[](2)) * Ex.operator[](2);

		if(Dot < CloseToZero)
		{
			CloseToZero = Dot;
			Q = i;
		}

    }

	Ex.operator/=(Ex.length());

	Vec3d S;
	double Dot = (basisPts[Q].operator[](0) - basisPts[0].operator[](0)) * Ex.operator[](0);
	Dot +=(basisPts[Q].operator[](1) - basisPts[0].operator[](1)) * Ex.operator[](1);
	Dot +=(basisPts[Q].operator[](2) - basisPts[0].operator[](2)) * Ex.operator[](2);

	S = Ex;
	S.operator*=(Dot);

	Vec3d T = basisPts[Q] - basisPts[0] - S;
	Vec3d Ey = T;
	Ey.operator/=(Ey.length());

	Vec3d R = basisPts[0];
	basisPts.clear();

	double uMin = 999999, uMax = -99999, vMin = 999999, vMax = -99999;
	vector<double> Us, Vs;

	for (int i=0; i < numPoints; i++) {

		double ui = (points[i].X - points[0].X) * Ex.operator[](0);
		ui += (points[i].Y - points[0].Y) * Ex.operator[](1);
		ui += (points[i].Z - points[0].Z) * Ex.operator[](2);

		if( ui < uMin)
			uMin = ui;
		if(ui > uMax)
			uMax = ui;

		double vi = (points[i].X - points[0].X) * Ey.operator[](0);
		vi += (points[i].Y - points[0].Y) * Ey.operator[](1);
		vi += (points[i].Z - points[0].Z) * Ey.operator[](2);

		if( vi < vMin)
			vMin = vi;
		if( vi > vMax)
			vMax = vi;

		Us.push_back(ui);
		Vs.push_back(vi);

    }

	uScale = uMax - uMin;
	vScale = vMax - vMin;

	for (int i=0; i < numPoints; i++) {

		Vec3d tmp = Vec3d((Us[i] - uMin)/uScale, (Vs[i] - vMin)/vScale, points[i].W);
		basisPts.push_back(tmp);

	}

    /******** BEGIN TODO ********/
printf("TODO: %s:%d\n", __FILE__, __LINE__); 

    /******** END TODO ********/
}



//
// TODO 3: ComputeHomography()
//		Computes the homography H from the plane specified by "points" to the image plane,
//		and its inverse Hinv.
//		If the plane is the reference plane (isRefPlane == true), don't convert the
//		coordinate system to the plane. Only do this for polygon patches where
//		texture mapping is necessary.
//		Coordinate system conversion is to be implemented in a separate routine
//		ConvertToPlaneCoordinate.
//		For more detailed explaination, see
//		http://www.cs.cornell.edu/courses/cs4670/2012fa/projects/p4/homography.pdf.
//
void ComputeHomography(CTransform3x3 &H, CTransform3x3 &Hinv, const vector<SVMPoint> &points, vector<Vec3d> &basisPts, bool isRefPlane)
{
    int i;
    int numPoints = (int) points.size();
    assert( numPoints >= 4 );

    basisPts.clear();
    if (isRefPlane) // reference plane
    {
        for (i=0; i < numPoints; i++) {
            Vec3d tmp = Vec3d(points[i].X, points[i].Y, points[i].W); // was Z, not W
            basisPts.push_back(tmp);
        }
    } 
    else // arbitrary polygon
    {
        double uScale, vScale; // unused in this function
        ConvertToPlaneCoordinate(points, basisPts, uScale, vScale);
    }

    // A: 2n x 9 matrix where n is the number of points on the plane
    //    as discussed in lecture
    int numRows = 2 * numPoints;
    const int numCols = 9;

    typedef Matrix<double, Dynamic, 9, RowMajor> MatrixType;
    MatrixType A = MatrixType::Zero(numRows, numCols);

    /******** BEGIN TODO ********/
    /* Fill in the A matrix for the call to MinEig */

	 for (int i = 0; i < numPoints; i++) {

		 double bx = double(points[i].u);
		 double by = double(points[i].v);
		 double ax = double(basisPts[i].operator[](0));
		 double ay = double(basisPts[i].operator[](1));

		 A(i * 2,0) = ax;     A(i * 2,1) = ay;     A(i * 2,2) = 1;       
		 A(i * 2,3) = 0;      A(i * 2,4) = 0;      A(i * 2,5) = 0;   
		 A(i * 2,6) = -bx*ax; A(i * 2,7) = -bx*ay; A(i * 2,8) = -bx;

 		 A(i * 2 + 1,0) = 0;      A(i * 2 + 1,1) = 0;      A(i * 2 + 1,2) = 0;   
		 A(i * 2 + 1,3) = ax;     A(i * 2 + 1,4) = ay;     A(i * 2 + 1,5) = 1;   
		 A(i * 2 + 1,6) = -by*ax; A(i * 2 + 1,7) = -by*ay; A(i * 2 + 1,8) = -by;

	 }

printf("TODO: %s:%d\n", __FILE__, __LINE__); 


    double eval, h[9];
    MinEig(A, eval, h);

    H[0][0] = h[0];
    H[0][1] = h[1];
    H[0][2] = h[2];

    H[1][0] = h[3];
    H[1][1] = h[4];
    H[1][2] = h[5];

    H[2][0] = h[6];
    H[2][1] = h[7];
    H[2][2] = h[8];

    /******** END TODO ********/

    // compute inverse of H
    if (H.Determinant() == 0)
        fl_alert("Computed homography matrix is uninvertible \n");
    else
        Hinv = H.Inverse();

    int ii;
    printf("\nH=[\n");
    for (ii=0; ii<3; ii++)
        printf("%e\t%e\t%e;\n", H[ii][0]/H[2][2], H[ii][1]/H[2][2], H[ii][2]/H[2][2]);
    printf("]\nHinv=[\n");

    for (ii=0; ii<3; ii++)
        printf("%e\t%e\t%e;\n", Hinv[ii][0]/Hinv[2][2], Hinv[ii][1]/Hinv[2][2], Hinv[ii][2]/Hinv[2][2]);

    printf("]\n\n");
}

