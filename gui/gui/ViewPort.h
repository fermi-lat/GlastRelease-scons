//     $Id$
//  Author: Toby Burnett
//


#ifndef VIEWPORT_H
#define VIEWPORT_H


#include "gui/GraphicsVector.h"

namespace gui {

class Draw2D;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class ViewPort {
    // Virtual base class for rendering 3-dimensional lines and points, with a
    // projection definition to a screen.
    // Itself implements a complete model that is consistent with OpenGL.
    //
    
public:
    ViewPort();
    ViewPort(const GraphicsVector*, float scale, const float* angles, const char* t=0);
    // default constructor, convenient for making copies
    
    
    virtual ViewPort& copy(const GraphicsVector* ref,
	float sc,
	const float * eulerAngles,
	const char* title=0);
    // create a new copy and set the parameters--subclass must implement to make a copy of itself
    
    virtual void set_canvas(Draw2D* new_canvas);
    // allow change of the canvas, if used by subclass
    
    virtual Draw2D* canvas()const;
    // returns pointer to current canvas object, if used by subclass
    
    // -----------------------------------
    // functions that affect the viewpoint
    // -----------------------------------
    virtual void reset();
    // restore defaults, initial orientation, anything subclass neeeds reset
    
    void rotateX(float);
    void rotateY(float);
    void rotateZ(float);
    void rotateEuler(const float*);
    // rotation
    
    
    void zoom(float z);
    void setZoom(float z);
    // adjust magnification, or zoom factor
    
    void moveIn(float = 0.05f);
    void moveOut(float = 0.05f);
    void setPerspective(float p);
    // adjust perspective, or view distance
    
    void setAspectRatio(float a);
    // set the height/width aspect ratio
    
    void setScale(float s);
    
    void setDetail(float);
    float detail()const;
    // set, report current detail factor
    
    
    // ---------------------------------------------------
    //   functions applying world<->window transformations
    // ---------------------------------------------------
    void transformToWindow(const GraphicsVector&, float&, float&);
    // project from world to window screen:
    
    GraphicsVector transformToWorld(float, float);
    // transform from window to world
    
    //  ----------------------------------------
    //  functions for setting drawing attributes
    //  ----------------------------------------
    virtual void set_rgb( float red, float green, float blue );
    virtual void set_line_width( float lwid );
    virtual void set_line_style( int s );
    
    virtual void set_col_index(int index);
    virtual void set_cap_style(int style);
    virtual void setJoinStyle(int style);
    virtual void setBkColIndex(int index);
    virtual void set_marker_size(float size);
    virtual void set_enhanced(int);
    // set the selected attribute on or off
    
    // -------------
    // draw commands
    // -------------
    virtual void set_quad(int quadrant, int highlighted);
    // setup the quadrant to draw into. 
    
    virtual void drawPL(const GraphicsVector*,int)=0;
    // draw a polyline 
    
    virtual void draw_markers(const GraphicsVector*,int)=0;
    // draw markers (points) at the specified points
    
    virtual void drawText(const char *text, const GraphicsVector& pnt, int ht);
    virtual void drawText(const char* text);
    // draw (horizontal) text at the given, or current point
    
    virtual void flush();
    // flush accumulated drawing commands. Called after all drawing commands
    
    // ---------------
    // VRML interface
    // --------------
    virtual void beginSeparator();
    virtual void endSeparator();
    // useful for VRML color, texture, etc.
    
    virtual void coordinate3(const GraphicsVector *points, int n)=0;
    // set the points to display
    
    virtual void indexedLineSet(const int* a, int n)=0;
    virtual void indexedFaceSet(const int* a, int n);
    // create polylines or set of faces from coordinate3 list

    // access
    float zoom()const;
    const char * title()const;
    void detail(float d); // set the detail level
    
 private:
     const char * _title;	// optional title
     
     float 
	 zoomFactor,  	// multiplication factor
	 aspectRatio, 	// screen aspectRatio
	 perspective, 	// inverse distance: 0 is infinity, 1 closest
	 scale;     	// scale factor to transform from world to model coordinates 
     const GraphicsVector* center; // center
     const float* angles; // euler angles expressing initial view orientation
     float R[16];     	// rotation matrix constructed from euler angles, subsequent rotations
     
     float _detail;
     static float defaultDetail;
};

//----------------------------------------------------------------------------
//            inlines
inline  void ViewPort::zoom(float z){zoomFactor*=z;}
inline float ViewPort::zoom()const{return zoomFactor;}
inline  void ViewPort::setZoom(float z){zoomFactor =z;}
inline  void ViewPort::moveIn(float f){perspective+=f;if(perspective>1)perspective=1;}
inline  void ViewPort::moveOut(float f){perspective-=f;if(perspective<0)perspective=0;}
inline  void ViewPort::setPerspective(float p){perspective=p;}
inline  void ViewPort::setAspectRatio(float a){aspectRatio=a;}
inline  void ViewPort::setScale(float s){scale= s;}           
inline void ViewPort::set_canvas(Draw2D*){}
inline Draw2D* ViewPort::canvas()const{return 0;}
inline void ViewPort::set_rgb( float , float , float  ){}
inline void ViewPort::set_line_width( float ){}
inline void ViewPort::set_line_style( int  ){}
inline void ViewPort::set_col_index(int){}
inline void ViewPort::set_cap_style(int ){}
inline void ViewPort::setJoinStyle(int){}
inline void ViewPort::setBkColIndex(int ){}
inline void ViewPort::set_marker_size(float ){}
inline void ViewPort::set_enhanced(int){}
inline void ViewPort::set_quad(int , int ){}
inline void ViewPort::drawText(const char *, const GraphicsVector& , int ){}
inline void ViewPort::drawText(const char*){}
inline void ViewPort::flush(){}
inline void ViewPort::beginSeparator(){}
inline void ViewPort::endSeparator(){}
inline void ViewPort::indexedFaceSet(const int* a, int n){indexedLineSet(a,n);}
inline const char* ViewPort::title()const{return _title;}
inline void ViewPort::detail(float d){_detail=d;}

} //namespace gui
#endif  //VIEWPORT_H


