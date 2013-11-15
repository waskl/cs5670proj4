/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #2:
 * ImgView.inl (included from ImgView.cpp)
 *		contains routines for computing the 3D position of points
 ***************************************************************/

//
// TODO 4: sameXY()
//		Computes the 3D position of newPoint using knownPoint
//		that has the same X and Y coordinate, i.e. is directly
//		below or above newPoint.
//		See lecture slide on measuring heights.
//
// HINT1: make sure to dehomogenize points when necessary
// HINT2: there is a degeneracy that you should look out for involving points already in line with the reference
// HINT3: make sure to get the sign of the result right, i.e. whether it is above or below ground
void ImgView::sameXY()
{
	if (pntSelStack.size() < 2)
	{
		fl_alert("Not enough points on the stack.");
		return;
	}

	SVMPoint &newPoint = *pntSelStack[pntSelStack.size() - 1];
	SVMPoint &knownPoint = *pntSelStack[pntSelStack.size() - 2];

	if( !knownPoint.known() )
	{
		fl_alert("Can't compute relative values for unknown point.");
		return;
	}

	if( refPointOffPlane == NULL )
	{
		fl_alert("Need to specify the reference height first.");
		return;
	}

	/******** BEGIN TODO ********/

	// See the lecture note on measuring heights
	// using a known point directly below the new point.
	 
	// printf("sameXY() to be implemented!\n");

	SVMPoint groundPoint(knownPoint); // Project known point onto ground for third Point
	groundPoint.Z = 0;
	groundPoint.u  = H[0][0] * groundPoint.X + H[0][1] * groundPoint.Y + H[0][2] * groundPoint.W;
	groundPoint.u /= H[2][0] * groundPoint.X + H[2][1] * groundPoint.Y + H[2][2] * groundPoint.W;

	groundPoint.v  = H[1][0] * groundPoint.X + H[1][1] * groundPoint.Y + H[1][2] * groundPoint.W;
	groundPoint.v /= H[2][0] * groundPoint.X + H[2][1] * groundPoint.Y + H[2][2] * groundPoint.W;

	groundPoint.w = 1;

	

	if(knownPoint.Z != 0)
	{

		double ty = (newPoint.v - groundPoint.v) / (knownPoint.v - groundPoint.v);
		double tx = (newPoint.u - groundPoint.u) / (knownPoint.u - groundPoint.u);
		double t = (tx + ty)/2;

		double t_k_image = sqrt( pow((newPoint.u - knownPoint.u),2) + pow((newPoint.v - knownPoint.v),2) ); // t - r
		double t_b_image = sqrt( pow((newPoint.u - groundPoint.u),2) + pow((newPoint.v - groundPoint.v),2) ); // t - b
		double vz_k_image = sqrt( pow((knownPoint.u - zVanish.u),2) + pow((knownPoint.v - zVanish.v),2) ); //vz - r
		double k_b_image = sqrt( pow((knownPoint.u - groundPoint.u),2) + pow((knownPoint.v - groundPoint.v),2) ); // r - b
		double vz_t_image = sqrt( pow((newPoint.u - zVanish.u),2) + pow((newPoint.v - zVanish.v),2) ); //vz - t
		double vz_b_image = sqrt( pow((groundPoint.u - zVanish.u),2) + pow((groundPoint.v - zVanish.v),2) ); //vz - b

		double Img_cross_ratio;

		if( t >= 0 && t <= 1) // point lies between R and B
		{
			printf("Point is between Reference and ground plane\n");
			Img_cross_ratio = (t_b_image * vz_k_image)/(k_b_image * vz_t_image); // == H/R
			Img_cross_ratio = 1 / Img_cross_ratio;

		}

		else if( t < 0 ) // point lies below B
		{
			printf("Point is below ground plane\n");
			Img_cross_ratio = (t_k_image * vz_b_image)/(t_b_image * vz_k_image); // == H/R
		}

		else if ( t > 0 ) // point lies above R
		{
			printf("Point is above Reference\n");
			Img_cross_ratio = (t_b_image * vz_k_image)/(k_b_image * vz_t_image); // == H/R
		}
	

		newPoint.X = knownPoint.X;
		newPoint.Y = knownPoint.Y;
		newPoint.Z = Img_cross_ratio * knownPoint.Z;

	}
	
	else if( knownPoint.X == refPointOffPlane->X && knownPoint.Y == refPointOffPlane->Y )
	{

		SVMPoint refPointOnGround;
		refPointOnGround.X = refPointOffPlane->X;
		refPointOnGround.Y = refPointOffPlane->Y;
		refPointOnGround.W = refPointOffPlane->W;
		refPointOnGround.Z = 0;

		refPointOnGround.u  = H[0][0] * refPointOnGround.X + H[0][1] * refPointOnGround.Y + H[0][2] * refPointOnGround.W;
		refPointOnGround.u /= H[2][0] * refPointOnGround.X + H[2][1] * refPointOnGround.Y + H[2][2] * refPointOnGround.W;

		refPointOnGround.v  = H[1][0] * refPointOnGround.X + H[1][1] * refPointOnGround.Y + H[1][2] * refPointOnGround.W;
		refPointOnGround.v /= H[2][0] * refPointOnGround.X + H[2][1] * refPointOnGround.Y + H[2][2] * refPointOnGround.W;

		refPointOnGround.w = 1;

		double ty = (newPoint.v - refPointOnGround.v) / (knownPoint.v - refPointOnGround.v);
		double tx = (newPoint.u - refPointOnGround.u) / (knownPoint.u - refPointOnGround.u);
		double t = (tx + ty)/2;

		double t_b_image = sqrt( pow(refPointOffPlane->u - refPointOnGround.u, 2) + pow(refPointOffPlane->v - refPointOnGround.v, 2) );
		double vz_b_image = sqrt( pow(zVanish.u - newPoint.u, 2) + pow(zVanish.v - newPoint.v, 2) );
		double unk_b_image = sqrt( pow(refPointOnGround.u - newPoint.u, 2) + pow(refPointOnGround.v - newPoint.v, 2) );
		double vz_t_image = sqrt( pow(zVanish.u - refPointOffPlane->u, 2) + pow(zVanish.v - refPointOffPlane->v, 2) );

		double ImgCrossRatio = (t_b_image * vz_b_image) / (unk_b_image * vz_t_image);

		newPoint.X = knownPoint.X;
		newPoint.Y = knownPoint.Y;
		newPoint.Z = refPointOffPlane->Z / ImgCrossRatio;

		if( t < 0 )
			newPoint.Z *= -1;

	}

	else
	{

		Vec3d VxVy;
		VxVy = cross( Vec3d(xVanish.u, xVanish.v, xVanish.w), Vec3d(yVanish.u, yVanish.v, yVanish.w) );
		
		SVMPoint refPointOnGround;
		refPointOnGround.X = refPointOffPlane->X;
		refPointOnGround.Y = refPointOffPlane->Y;
		refPointOnGround.W = refPointOffPlane->W;
		refPointOnGround.Z = 0;

		refPointOnGround.u  = H[0][0] * refPointOnGround.X + H[0][1] * refPointOnGround.Y + H[0][2] * refPointOnGround.W;
		refPointOnGround.u /= H[2][0] * refPointOnGround.X + H[2][1] * refPointOnGround.Y + H[2][2] * refPointOnGround.W;

		refPointOnGround.v  = H[1][0] * refPointOnGround.X + H[1][1] * refPointOnGround.Y + H[1][2] * refPointOnGround.W;
		refPointOnGround.v /= H[2][0] * refPointOnGround.X + H[2][1] * refPointOnGround.Y + H[2][2] * refPointOnGround.W;

		refPointOnGround.w = 1;

		Vec3d RcrossRefPtGrnd;
		RcrossRefPtGrnd = cross( Vec3d(knownPoint.u, knownPoint.v, knownPoint.w), Vec3d(refPointOnGround.u, refPointOnGround.v, refPointOnGround.w));

		Vec3d V = cross ( VxVy, RcrossRefPtGrnd );
		Vec3d UnknownV = cross( Vec3d( newPoint.u, newPoint.v, newPoint.w), V );
		Vec3d RefHtLn = cross( Vec3d(refPointOffPlane->u,refPointOffPlane->v,refPointOffPlane->w), Vec3d(refPointOnGround.u, refPointOnGround.v, refPointOnGround.W));
		Vec3d IntMdPt = cross( RefHtLn, UnknownV );

		IntMdPt[0] = IntMdPt[0] / IntMdPt[2];
		IntMdPt[1] = IntMdPt[1] / IntMdPt[2];
		IntMdPt[2] = IntMdPt[2] / IntMdPt[2];

		double ty = (IntMdPt[1] - refPointOnGround.v) / (refPointOffPlane->v - refPointOnGround.v);
		double tx = (IntMdPt[0] - refPointOnGround.u) / (refPointOffPlane->u - refPointOnGround.u);
		double t = (tx + ty)/2;

		double t_b_image = sqrt( pow(refPointOffPlane->u - refPointOnGround.u, 2) + pow(refPointOffPlane->v - refPointOnGround.v, 2) );
		double vz_b_image = sqrt( pow(zVanish.u - IntMdPt[0], 2) + pow(zVanish.v - IntMdPt[1], 2) );
		double unk_b_image = sqrt( pow(refPointOnGround.u - IntMdPt[0], 2) + pow(refPointOnGround.v - IntMdPt[1], 2) );
		double vz_t_image = sqrt( pow(zVanish.u - refPointOffPlane->u, 2) + pow(zVanish.v - refPointOffPlane->v, 2) );

		double ImgCrossRatio = (t_b_image * vz_b_image) / (unk_b_image * vz_t_image);

		newPoint.X = knownPoint.X;
		newPoint.Y = knownPoint.Y;
		newPoint.Z = refPointOffPlane->Z / ImgCrossRatio;

		if( t < 0 )
			newPoint.Z *= -1;

	}
	
	
printf("TODO: %s:%d\n", __FILE__, __LINE__); 

	/******** END TODO ********/

	newPoint.known(true);

	printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );

	redraw();
}



//
// TODO 5: sameZPlane()
//		Compute the 3D position of newPoint using knownPoint
//		that lies on the same plane and whose 3D position is known.
//		See the man on the box lecture slide.
//		If newPoint is on the reference plane (Z==0), use homography (this->H, or simply H) directly.
//
// HINT: For this function, you will only need to use the three vanishing points and the reference homography 
//       (in addition to the known 3D location of knownPoint, and the 2D location of newPoint)
void ImgView::sameZPlane()
{
	if (pntSelStack.size() < 2)
	{
		fl_alert("Not enough points on the stack.");
		return;
	}

	SVMPoint &newPoint = *pntSelStack[pntSelStack.size() - 1];
	SVMPoint &knownPoint = *pntSelStack[pntSelStack.size() - 2];

	if( !knownPoint.known() )
	{
		fl_alert("Can't compute relative values for unknown point.");
		return;
	}

	/******** BEGIN TODO ********/
	
	Vec3d b0 = Vec3d(newPoint.u, newPoint.v, newPoint.w);

	if(knownPoint.Z != 0)
	{

		Vec3d Vx = Vec3d(xVanish.u, xVanish.v, xVanish.w);
        Vec3d Vy = Vec3d(yVanish.u, yVanish.v, yVanish.w);
        Vec3d Vz = Vec3d(zVanish.u, zVanish.v, zVanish.w);
		
		Vec3d newPt = Vec3d(newPoint.u, newPoint.v, newPoint.w);
		Vec3d knownPt = Vec3d(knownPoint.u, knownPoint.v, knownPoint.w);
		
		Vec3d newPtknownPtLine = cross (newPt, knownPt);
		Vec3d VxVy = cross( Vx, Vy );

		Vec3d VxVyIntersection = cross( VxVy, newPtknownPtLine );

		SVMPoint knownPtGround(knownPoint);
		knownPtGround.W = 1;
		knownPtGround.Z = 0;

		knownPtGround.u  = H[0][0] * knownPtGround.X + H[0][1] * knownPtGround.Y + H[0][2] * knownPtGround.W;
		knownPtGround.u /= H[2][0] * knownPtGround.X + H[2][1] * knownPtGround.Y + H[2][2] * knownPtGround.W;

		knownPtGround.v  = H[1][0] * knownPtGround.X + H[1][1] * knownPtGround.Y + H[1][2] * knownPtGround.W;
		knownPtGround.v /= H[2][0] * knownPtGround.X + H[2][1] * knownPtGround.Y + H[2][2] * knownPtGround.W;

		knownPtGround.w = 1;

		Vec3d b1(knownPtGround.u, knownPtGround.v, 1.0);

		Vec3d b1VLine = cross (b1, VxVyIntersection);

		Vec3d VznewPt = cross(newPt, Vz);

		b0 = cross( VznewPt, b1VLine);
		b0[0] = b0[0]/b0[2];
		b0[1] = b0[1]/b0[2];
		b0[2] = b0[2]/b0[2];


	}
	
	double newX, newY;
	newX  = Hinv[0][0] * b0[0] + Hinv[0][1] * b0[1] + Hinv[0][2] * b0[2];
	newX /= Hinv[2][0] * b0[0] + Hinv[2][1] * b0[1] + Hinv[2][2] * b0[2];
																		  
	newY  = Hinv[1][0] * b0[0] + Hinv[1][1] * b0[1] + Hinv[1][2] * b0[2];
	newY /= Hinv[2][0] * b0[0] + Hinv[2][1] * b0[1] + Hinv[2][2] * b0[2];

	newPoint.X = newX;
	newPoint.Y = newY;
	newPoint.Z = knownPoint.Z;
	newPoint.W = 1.0;

printf("TODO: %s:%d\n", __FILE__, __LINE__); 

	/******** END TODO ********/

	newPoint.known(true);

	printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );

	redraw();
}


