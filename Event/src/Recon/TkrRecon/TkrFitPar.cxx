//#define TIMETRIAL 1
#undef TIMETRIAL

#ifdef TIMETRIAL
#include "TkrFitPar.h"
#else /* Standard Gleam */
#include "Event/Recon/TkrRecon/TkrFitPar.h"
#include "geometry/Ray.h"
#endif /* TIMETRIAL */
#include "CLHEP/Matrix/Vector.h"

using namespace Event;

TkrFitPar::TkrFitPar () 
{
  m_xpos = 0.0;
  m_xslope = 0.0;
  m_ypos = 0.0;
  m_yslope = 0.0;

}

TkrFitPar::TkrFitPar(const TkrFitPar &p)
{
  m_xpos = p.m_xpos;
  m_xslope = p.m_xslope;
  m_ypos = p.m_ypos;
  m_yslope = p.m_yslope;

}


TkrFitPar::TkrFitPar(const HepVector &p)
{
  m_xpos = p(1);
  m_xslope = p(2);
  m_ypos = p(3);
  m_yslope = p(4);

}

#ifndef TIMETRIAL
TkrFitPar::TkrFitPar(const Ray &ray) 
{

  Vector dir = ray.direction();
  Point  x0  = ray.position();
  double x_slope = dir.x()/dir.z();
  double y_slope = dir.y()/dir.z(); 
  double x_inter = x0.x();
  double y_inter = x0.y();
  m_xpos = x_inter;
  m_xslope = x_slope;
  m_ypos = y_inter;
  m_yslope = y_slope;

}

TkrFitPar::TkrFitPar(const Ray *ray) 
{
  Vector dir = ray->direction();
  Point  x0  = ray->position();
  double x_slope = dir.x()/dir.z();
  double y_slope = dir.y()/dir.z(); 
  double x_inter = x0.x();
  double y_inter = x0.y();
  m_xpos = x_inter;
  m_xslope = x_slope;
  m_ypos = y_inter;
  m_yslope = y_slope;
}

#endif TIMETRIAL


// Create with explicit values
TkrFitPar::TkrFitPar(double ax, double sx, double ay,double sy) 
{ 
  m_xpos = ax; 
  m_xslope = sx; 
  m_ypos = ay; 
  m_yslope = sy;

}

// Access methods for individual fit parameters
double TkrFitPar::getXPosition() const {return m_xpos;}
double TkrFitPar::getXSlope()    const {return m_xslope;}
    
double TkrFitPar::getYPosition() const {return m_ypos;}
double TkrFitPar::getYSlope()    const {return m_yslope;}  
    

//operator overload for element access TkrFitPar(element)
//!NOTE index starts with (1). In case of index error returns (1)
const double & TkrFitPar::operator() (const int &element) const
{
  if (element == 1 ) {
    return m_xpos;
  } else if (element == 2) {
    return m_xslope;
  } else if (element == 3) {
    return m_ypos;
  } else if (element == 4) {
    return m_yslope;
  } else {
    return m_xpos;
  }  
      
      
}

//!NOTE index starts with (1). In case of index error returns (1)
double & TkrFitPar::operator() (const int &element)
{
  if (element == 1 ) {
    return m_xpos;
  } else if (element == 2) {
    return m_xslope;
  } else if (element == 3) {
    return m_ypos;
  } else if (element == 4) {
    return m_yslope;
  } else {
    return m_xpos;
  }  
 
}

//operator overload for element access TkrFitPar(element)
//!NOTE index starts with [0]. In case of index error returns [0]
const double & TkrFitPar::operator[] (const int &element) const
{
  if (element == 0 ) {
    return m_xpos;
  } else if (element == 1) {
    return m_xslope;
  } else if (element == 2) {
    return m_ypos;
  } else if (element == 3) {
    return m_yslope;
  } else {
    return m_xpos;
  }  
      
      
}

//!NOTE index starts with [0]. In case of index error returns [0]
double & TkrFitPar::operator[] (const int &element)
{
  if (element == 0 ) {
    return m_xpos;
  } else if (element == 1) {
    return m_xslope;
  } else if (element == 2) {
    return m_ypos;
  } else if (element == 3) {
    return m_yslope;
  } else {
    return m_xpos;
  }  
 
}

//operator overload for TkrFitPar+TkrFitPar
const TkrFitPar TkrFitPar::operator +(const TkrFitPar& A) const
{
    
  double xpos = m_xpos + A.m_xpos;
  double xslope = m_xslope + A.m_xslope;
  double ypos = m_ypos + A.m_ypos;
  double yslope = m_yslope + A.m_yslope;

  return TkrFitPar(xpos,xslope,ypos,yslope);

}

//operator overload for TkrFitPar-TkrFitPar
const TkrFitPar TkrFitPar::operator -(const TkrFitPar& A) const
{
    
  double xpos = m_xpos - A.m_xpos;
  double xslope = m_xslope - A.m_xslope;
  double ypos = m_ypos - A.m_ypos;
  double yslope = m_yslope - A.m_yslope;

  return TkrFitPar(xpos,xslope,ypos,yslope);

}

//operator overload for TkrFitPar*TkrFitPar
const double TkrFitPar::operator *(const TkrFitPar& A) const
{
    
  return m_xpos * A.m_xpos + m_xslope * A.m_xslope
    + m_ypos * A.m_ypos + m_yslope * A.m_yslope;


}


