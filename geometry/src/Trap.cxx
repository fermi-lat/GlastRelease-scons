//  $Id$
//

#include "geometry/Trap.h"

static Vector N(const Vector& A, const Vector& B, const Vector& C)
{
   // define a plane with vectors to three (non-colinear) Vectors
   Vector AB = A-B,
          AC = A-C,
   		  n  = (AB.cross(AC)).unit();
   return (n*A) * n;
}

Trap::Trap( double h, double th,
	    double w1, double bl1, double tl1, double th1,
	    double w2, double bl2, double tl2 )
: Hexahedron()
{
//  Make a Trap with the following 9 parameters:
//	 1) the full length in z (h)
//	 2) the polar angle of the center of the -z face to the center of
//		the +z face (th) (-PI/2 < th < +PI/2)
//	 3) the full length in y of the -z face (w1)
//	 4) the full length in x of the low y edge of the -z face (bl1)
//	 5) the full length in x of the high y edge of the -z face (tl1)
//	 6) the angle with respect to the y-axis from the center of the
//		low y edge to the center of the high y edge of the -z face
//		(th1) (-PI/2 < th1 < +PI/2)
//	 7) the full length in y of the +z face (w2)
//	 8) the full length in x of the low y edge of the +z face (bl2)
//	 9) the full length in x of the high y edge of the +z face (tl2)
//  The center of the Trap is at the origin of the global coordinate system
//  and the axes are lined up with the global x (length), y (width), and
//  z (height) axes
	if(   ( h <= 0.0 )
	   || ( w1 <= 0.0 ) || ( bl1 <= 0.0 ) || ( tl1 <= 0.0 )
	   || ( w2 <= 0.0 ) || ( bl2 <= 0.0 ) || ( tl2 <= 0.0 )
	   || fabs( th ) >= 0.5 * M_PI
	   || fabs( th1) >= 0.5 * M_PI
	  )	FATAL("Illegal parameters for a Trap");

	double hh = 0.5 * h;
	double hw1 = 0.5 * w1;
	double hw2 = 0.5 * w2;
	double hbl1 = 0.5 * bl1;
	double hbl2 = 0.5 * bl2;
	double htl1 = 0.5 * tl1;
	double htl2 = 0.5 * tl2;
	double hw1t = hw1 * tan( th1 );
	double hht = hh * tan( th );
    //  Setup the vertices of the Trap
	Vector A (  hht + htl2 + hw1t,   hw2,  hh );
	Vector B (  hht - htl2 + hw1t,   hw2,  hh );
	Vector C (  hht + hbl2 - hw1t,  -hw2,  hh );
	Vector D (  hht - hbl2 - hw1t,  -hw2,  hh );
	Vector E ( -hht + htl1 + hw1t,   hw1, -hh );
	Vector F ( -hht - htl1 + hw1t,   hw1, -hh );
	Vector G ( -hht + hbl1 - hw1t,  -hw1, -hh );
	Vector H ( -hht - hbl1 - hw1t,  -hw1, -hh );
//  Construct six planes
    addPlane(N(A,C,G));	 //(+x)
	addPlane(N(D,B,H));	 //(-x)
	addPlane(N(B,A,F));	 //(+y)
	addPlane(N(C,D,G));	 //(-y)
	addPlane(N(A,B,C));	 //(+z)
	addPlane(N(F,E,H));	 //(-z)
}

	double Trap::height() const{return 0 ;}	  // put these off for now, or until Bill complains
	double Trap::theta() const{return 0 ;}
	double Trap::width1() const{return 0 ;}
	double Trap::lengthLow1() const{return 0 ;}
	double Trap::lengthHigh1() const {return 0 ;}
	double Trap::theta1() const {return 0 ;}
	double Trap::width2() const {return 0 ;}
	double Trap::lengthLow2() const {return 0 ;}
	double Trap::lengthHigh2() const {return 0 ;}

void Trap::printOn( std::ostream& os ) const
{
     Volume::printOn(os);
	 os
	   << "height "       << height()
	   << ", theta "      << theta()
	   << "\nwidth_1 "    << width1()
	   << ", length_low_1 " << lengthLow1()
	   << ", length_high_1 " << lengthHigh1()
	   << ", theta_1 "    << theta1()
	   << "\nwidth_2 "    << width2()
	   << ", length_low_2 " << lengthLow2()
	   << ", length_high_2 " << lengthHigh2()
	   << "\n" ;
}


