//     $Id$
//  Author: Toby Burnett
//
//

#ifdef _MSC_VER
#pragma warning(disable: 4244)  // conversion of float to double warning
#endif

#include "gui/ViewPort.h"
#include "gui/GraphicsVector.h"
#include <cmath>

namespace gui {
// these for convenience of subclasses that do no need to project (like VRML)
float ViewPort::defaultDetail= 0.01f;

static GraphicsVector currentPoint;
static GraphicsVector defaultOrigin(0,0,0);
static float defaultEulers[3]={0,0,0};
static const char* defaultTitle="";

ViewPort::ViewPort(const GraphicsVector* _center, float sc,const float* _angles, const char* t)
: _title(t)
, aspectRatio(1)
, scale(sc)
, center(_center)
, angles(_angles)
{
  reset();
}

// default constructor for subclasses
ViewPort::ViewPort()
: _title(defaultTitle)
, aspectRatio(1)
, scale(100.)
, center(&defaultOrigin)
, angles(defaultEulers)
{
    reset();
}
void ViewPort::reset()  
{ 
   _detail = defaultDetail;
   zoomFactor =1; perspective=0;
 
   float identity[]={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
   for(int i=0; i<16; i++){ R[i] = identity[i];}
   rotateEuler(angles);

   perspective = angles[3];

}

void
ViewPort::transformToWindow(const GraphicsVector& pt, float& winx,float& winy)
{
	float deltax = pt.x() - center->x();   //ugly, but efficient, I hope
	float deltay = pt.y() - center->y();
	float deltaz = pt.z() - center->z();
	float t = zoomFactor/scale;
	winx = (R[0]*deltax + R[4]*deltay + R[8]*deltaz)*t;
	winy = (R[1]*deltax + R[5]*deltay + R[9]*deltaz)*t*aspectRatio;

	if( perspective > 0.){
	   float winz = (R[2]*deltax + R[6]*deltay + R[10]*deltaz)*t;
	   float tt = 1.- winz*perspective;
	   if( tt > 1.e-3){
		winx /= tt; winy /= tt;
	    }
	}
}
 
GraphicsVector 
ViewPort::transformToWorld(float winx, float winy)
{
   return GraphicsVector
          ( (winx*R[0] + winy*R[1]) * scale/zoomFactor
          , (winx*R[4] + winy*R[5]) * scale/zoomFactor/aspectRatio
	  , (winx*R[8] + winy*R[9]) * scale/zoomFactor 
	  ); 
} 


float ViewPort::detail()const{  return _detail;}
void ViewPort::setDetail(float f){_detail=f;}
//----------------------------------------------------------------------
//                      Rotations--very efficient

void
ViewPort::rotateEuler(const float * angles)
{
        rotateZ(angles[0]);
        rotateY(angles[1]);
        rotateZ(angles[2]);
}

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

static const double d2r = M_PI/180;
static double c, s;

static inline void
rotate2(float& a, float &b)
{
  float t = c*a-s*b;
  b = s*a+c*b;
  a =t;
}  
void
ViewPort::rotateX(float theta)
{
      if( theta == 0. ) return;
      c = cos(theta*=d2r);
      s = sin(theta);
      rotate2(R[1],R[2]);
      rotate2(R[5],R[6]);
      rotate2(R[9],R[10]);
}
void
ViewPort::rotateY(float theta)
{
   if( theta ==0. ) return;
      c = cos(theta*=d2r);
      s =-sin(theta);
      rotate2(R[0],R[2]);
      rotate2(R[4],R[6]);
      rotate2(R[8],R[10]);
}

void
ViewPort::rotateZ(float theta)
{    if( theta ==0. ) return;
      c = cos(theta*=d2r);
      s = sin(theta);
      rotate2(R[0],R[1]);
      rotate2(R[4],R[5]);
      rotate2(R[8],R[9]);
}    
   // dummy! Better to put this inside the next function, but KAI chokes.
    class ViewPortDummy : public ViewPort { public:
    virtual void drawPL(const GraphicsVector*,int){};
    virtual void draw_markers(const GraphicsVector*,int){};
    virtual void coordinate3(const GraphicsVector *points, int n){};
    virtual void indexedLineSet(const int* a, int n){};
    };
 
ViewPort& ViewPort::copy(const GraphicsVector*, 	float, const float *,const char*)
{
    return *new ViewPortDummy;
}


} // namespace gui
