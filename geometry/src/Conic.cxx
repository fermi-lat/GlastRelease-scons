//     $Id$
//  Author: T. Burnett
// Project: geometry
//
// Conic class implementation

#include "geometry/Conic.h"

inline static double sqr(double x){return x*x;}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Conic::Conic( const Point& origin, const Vector& axis,
             double radius, double slope )
//----------------------------------------------------
   : Surface(origin ,axis.unit())
, m_slope(slope)
, m_radius(radius)

{}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double
Conic::how_near( const Point& x ) const
//------------------------------------
//  Distance from the point x to the nearest point on the Conic.
//  The distance will be positive if the point is inside the Conic,
//  negative if the point is outside.
{
    Vector r = x - origin();
    double z = r*axis(),
       rho = sqrt( r.mag2() - sqr(z) ),
       c = sqrt(1+sqr(m_slope));

    return (radius()>0?  radius(z)-rho  : rho+radius(z))/c ;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Vector
Conic::normal( const Point& x )const
//----------------------------------
{
    Vector r = x-origin();
    Vector rhohat = (r - (r*axis())*axis()).unit(); // unit vector in rho direction
    double b = radius()>0? m_slope : -m_slope;		// actual slope

    Vector n = sqrt(1-sqr(b))*rhohat - b* axis();
    return radius()>0 ? n : -n;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double
Conic::distance(const Point& x, const Vector& vhat, int inout)const
//-----------------------------------------------------------------
//   solve: | rperp(s) | = R(rpar(s)),
//   where rperp(s) = r(s)-ahat*(r(s)*ahat)
//         rpar(s) = r(s)*ahat;
//         r(s) = x-origin() + vhat*s
//	    R(z) = R0 + alpha*z
{

    Vector  ahat = axis(),
	    r = x - origin();
    double  rpar = r*ahat,
	    vpar = vhat*ahat,
	    alpha= m_slope,
	    R = radius(rpar);
    Vector  rperp = r - rpar*ahat,
	    vperp = vhat-vpar*ahat;
	
    double a = vperp.mag2()-sqr(alpha*vpar),
	   b = vperp*rperp-alpha*vpar*R,
	   c = rperp.mag2() - sqr(R),
	   disc = b*b-a*c,
	   s;
    // now want appropriate root of a*s^2 + 2*b*s + c = 0

    if( disc<=0 ) return FLT_MAX;   // misses
    if( radius()<0 ) inout *=-1;    // reverse enter/leave sense
    int leaving = inout==1;
    if( b<0 )  {
	// velocity toward axis
	s = (leaving) ? (-b+sqrt(disc))/a :
                          c/(-b+sqrt(disc));
    } else {
	// heading away from axis: can only leave
	s = (leaving) ? c/(-b-sqrt(disc)) : FLT_MAX;
    }
    // finally see if the solution is invalid, i.e., past the apex of the cone
    if( m_slope!=0 && (rpar + s*vpar) < -m_radius/m_slope ){
        //  std::cerr << "s="<< s<< ", rpar(s)= " << rpar+s*vpar << ", apex at " <<- m_radius/m_slope<<std::endl ;;
       s = FLT_MAX;
    }
    return s ;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Conic::printOn( std::ostream& os ) const
//--------------------------------------
{
    	os << "Conic surface with origin: " << origin()
	   << ",\t radius: " << radius()
	   << ",\t axis: " << axis()
	   << ",\t slope: "<< m_slope << "\n";
}



