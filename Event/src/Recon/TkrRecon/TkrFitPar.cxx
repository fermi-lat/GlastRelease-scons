
// $Header$

//----------------------------------------------------------------------
//    
//    Implementation of the Kalman Filter Functions
//                TkrFitPar 
//
//      Original due to Jose Hernando-Angel circa 1997-1999
//      Re-written to combine both X and Y projections (2001) 
//      
//-----------------------------------------------------------------------

#include "GlastEvent/Recon/TkrRecon/TkrFitPar.h"
#include "geometry/Ray.h"

using namespace TkrRecon;

TkrFitPar::TkrFitPar(const Ray &ray) : HepVector(4)
{
    Vector dir = ray.direction();
    Point  x0  = ray.position();
    double x_slope = dir.x()/dir.z();
    double y_slope = dir.y()/dir.z(); 
    double x_inter = x0.x();
    double y_inter = x0.y();
    operator[](0) = x_inter;
    operator[](1) = x_slope;
    operator[](2) = y_inter;
    operator[](3) = y_slope;
}

TkrFitPar::TkrFitPar(const Ray *ray) : HepVector(4)
{
    Vector dir = ray->direction();
    Point  x0  = ray->position();
    double x_slope = dir.x()/dir.z();
    double y_slope = dir.y()/dir.z(); 
    double x_inter = x0.x();
    double y_inter = x0.y();
    operator[](0) = x_inter;
    operator[](1) = x_slope;
    operator[](2) = y_inter;
    operator[](3) = y_slope;
}

