/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #3:
 * ImgViewBox.cpp
 *		routines for finding the corners of an axis-aligned
 *      box (in 2D and in 3D)
 **************************************************************/
#pragma warning(disable : 4996)

#include "ImgView.h"

//
// TODO 6: solveForOppositeCorners()
//     Given the 2D positions of two corners of a rectangular face parallel to the XZ plane, compute
//     the 2D positions of the other two corners
void ImgView::solveForOppositeCorners(double u0, double v0, double u2, double v2,
                                      double &u1, double &v1, double &u3, double &v3)
{
    /* Vanishing points must be known */
    assert(xVanish.known() && yVanish.known() && zVanish.known());    


    /******** BEGIN TODO ********/ 
    // Given the 2D positions of corners p0 and p2 of the face, compute the 2D positions of p1 and p3
    // Remember that this face is on a plane perpendicular to the plane x=0
    // Store the results in variables 'u1, v1' and 'u3, v3'

//printf("TODO: %s:%d\n", __FILE__, __LINE__); 

	double slopeOne = (xVanish.v - v0) / (xVanish.u - u0);
	double slopeTwo = (zVanish.v - v2) / (zVanish.u - u2);

	u1 = ((slopeOne * xVanish.u) + zVanish.v - xVanish.v - (slopeTwo * zVanish.u)) / (slopeOne - slopeTwo);
	v1 = slopeOne * (u1 - xVanish.u) + xVanish.v;

	slopeOne = (zVanish.v - v0) / (zVanish.u - u0);
	slopeTwo = (xVanish.v - v2) / (xVanish.u - u2);

	u3 = ((slopeOne * zVanish.u) + xVanish.v - zVanish.v - (slopeTwo * xVanish.u)) / (slopeOne - slopeTwo);
	v3 = slopeOne * (u3 - zVanish.u) + zVanish.v;


    /********* END TODO ********/
}

//
// TODO 7: solveForOppositeFace()
//     Given the 2D positions of one rectangular face parallel to the XZ plane, 
//     compute the 2D positions of a parallel face being swept out from it.
//     The mouse position is given; one of the lines on the parallel face should pass
//     through the mouse position
void ImgView::solveForOppositeFace(SVMSweep *sweep, double imgX, double imgY,
                                   Vec3d &p4_out, Vec3d &p5_out, Vec3d &p6_out, Vec3d &p7_out)
{
    SVMPolygon *poly = sweep->poly;

    if (poly == NULL)
        return;

    // Get the four existing points
    SVMPoint *n0, *n1, *n2, *n3;
    poly->getFourPoints(&n0, &n1, &n2, &n3);

    Vec3d p0(n0->u, n0->v, n0->w);
    Vec3d p1(n1->u, n1->v, n1->w);
    Vec3d p2(n2->u, n2->v, n2->w);
    Vec3d p3(n3->u, n3->v, n3->w);

    Vec3d pMouse(imgX, imgY, 1.0);

    /******** BEGIN TODO ********/
    // Find the 2D image positions of box corners p4, p5, p6, p7, as described on the webpage.  
    // You will compute these positions using the known corners of the box (p0, p1, p2, p3, defined above)
    // and the vanishing points.
    // The line through points p4 and p5 will go through the mouse position, pMouse
    // Store the results in variables p4, p5, p6, and p7.
	Vec3d p4, p5, p6, p7;

//printf("TODO: %s:%d\n", __FILE__, __LINE__); 
	
	double slopeOne = (yVanish.v - p1[1]) / (yVanish.u - p1[0]);
	double slopeTwo = (xVanish.v - pMouse[1]) / (xVanish.u - pMouse[0]);

	p5[0] = ((slopeOne * yVanish.u) + xVanish.v - yVanish.v - (slopeTwo * xVanish.u)) / (slopeOne - slopeTwo);
	p5[1] = slopeOne * (p5[0] - yVanish.u) + yVanish.v;
	p5[2] = 1.0;

	slopeOne = (yVanish.v - p0[1]) / (yVanish.u - p0[0]);
	slopeTwo = (xVanish.v - p5[1]) / (xVanish.u - p5[0]);

	p4[0] = ((slopeOne * yVanish.u) + xVanish.v - yVanish.v - (slopeTwo * xVanish.u)) / (slopeOne - slopeTwo);
	p4[1] = slopeOne * (p4[0] - yVanish.u) + yVanish.v;
	p4[2] = 1.0;

	slopeOne = (zVanish.v - p4[1]) / (zVanish.u - p4[0]);
	slopeTwo = (yVanish.v - p3[1]) / (yVanish.u - p3[0]);

	p7[0] = ((slopeOne * zVanish.u) + yVanish.v - zVanish.v - (slopeTwo * yVanish.u)) / (slopeOne - slopeTwo);
	p7[1] = slopeOne * (p7[0] - zVanish.u) + zVanish.v;
	p7[2] = 1.0;

	slopeOne = (yVanish.v - p2[1]) / (yVanish.u - p2[0]);
	slopeTwo = (xVanish.v - p7[1]) / (xVanish.u - p7[0]);

	p6[0] = ((slopeOne * yVanish.u) + xVanish.v - yVanish.v - (slopeTwo * xVanish.u)) / (slopeOne - slopeTwo);
	p6[1] = slopeOne * (p6[0] - yVanish.u) + yVanish.v;
	p6[2] = 1.0;

    /******** END TODO ********/

    p4_out = p4;
    p5_out = p5;
    p6_out = p6;
    p7_out = p7;
}

//
// TODO 8: find3DPositionsBox()
//    Find the 3D positions of the 8 corners of the box.  The 3D position of points[0] is known.
void ImgView::find3DPositionsBox(SVMPoint *points[8]) 
{
    /******** BEGIN TODO ********/
    // Implement this function.  You will compute the 3D positions of the corners of the box
    // using their 2D positions and the known 3D position of points[0].  Store the results in 
    // points[1] through points[7].  You can and should use the sameXY and sameZ routines that 
	// you need to implement.  For that to work, you will need to push and pop points from
	// pntSelStack.  There are multiple ways to implement this function.

//printf("TODO: %s:%d\n", __FILE__, __LINE__); 
	pntSelStack.push_back(points[0]);
	pntSelStack.push_back(points[1]);
	sameXY();
	pntSelStack.pop_back();
	pntSelStack.push_back(points[4]);
	sameXY();
	pntSelStack.pop_back();
	pntSelStack.push_back(points[5]);
	sameXY();
	pntSelStack.pop_back();

	pntSelStack.push_back(points[3]);
	sameZPlane();
	pntSelStack.push_back(points[2]);
	sameXY();
	pntSelStack.pop_back();	
	pntSelStack.push_back(points[6]);
	sameXY();
	pntSelStack.pop_back();
	pntSelStack.push_back(points[7]);
	sameXY();
	pntSelStack.pop_back();		//pop 7
	pntSelStack.pop_back();		//pop 3
	pntSelStack.pop_back();		//pop 0
	
	/********* END TODO ********/
}


