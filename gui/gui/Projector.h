//     $Id$
//  Author: Toby Burnett
//
//
#ifndef PROJECTOR_H
#define PROJECTOR_H

#include "gui/ViewPort.h"
#include "gui/Draw2D.h"

namespace gui {

class GraphicsVector;

// Renders 3-dimensional lines and points, with a projection definition,
// to a screen
//
//  The projection is implemented by the superclass
//  actual drawing to a screen is implemented by the associated device-depended
//  Draw2D
class Projector : public ViewPort
{
public:
    Projector(Draw2D* d);
    // default constructor, convenient for making copies
    
    Projector(Draw2D*d, const GraphicsVector* ref,
	float sc,
	const float * eulerAngles,
	const char* title=0);
    // constructor: specify the 2-d drawer, the reference point,
    // the scales, and the initial orientation
    
    ViewPort& copy( const GraphicsVector* ref,
	float sc,
	const float * eulerAngles,
	const char* title=0);
    // create a new copy and set the parameters
    
    
    virtual void set_canvas(Draw2D* newDrawer);
    // allow change of the canvas
    
    virtual Draw2D* canvas()const;
    // returns pointer to current canvas object
    
    //  ----------------------------------------
    //  functions for setting drawing attributes
    //  ----------------------------------------
    void set_rgb( float red, float green, float blue );
    void set_line_width( float lwid );
    void set_line_style( int s );
    
    void set_col_index(int index);
    void set_cap_style(int style);
    void set_marker_size(float size);
    void set_enhanced(int);
    // set the selected attribute on or off
    
    // -------------
    // draw commands
    // -------------
    void set_quad(int quadrant, int highlighted);
    
    void drawPL(const GraphicsVector*,int);
    // draw a polyline 
    
    void drawText(const char *text, const GraphicsVector& pnt, int ht);
    void drawText(const char* text);
    // draw (horizontal) text at the given, or current point
    
    void draw_markers(const GraphicsVector*,int);
    // draw markers (points) at the specified points
    
    // ---------------
    // VRML interface
    // --------------
    void beginSeparator(){};
    void endSeparator(){};
    // useful for VRML color, texture, etc.
    
    void coordinate3(const GraphicsVector *points, int n);
    // set the points to display
    
    void indexedLineSet(const int* a, int n);
    void indexedFaceSet(const int* a, int n);
    // create polylines or set of faces from coordinate3 list
    
    
    
    void flush();
    // flush accumulated drawing commands
    
    
private:
    Draw2D* m_canvas;  	// pointer to the 2-d draw object
    
};

// ------------- inline functions ----------------------------------

inline void Projector::set_canvas(Draw2D* newDrawer){m_canvas=newDrawer;}
inline Draw2D* Projector::canvas()const{return m_canvas;}

inline void Projector::set_rgb( float red, float green, float blue ){
    m_canvas->set_rgb(red, green, blue);}

inline void Projector::set_line_width( float lwid ){
    m_canvas->set_line_width(lwid);}

inline void Projector::set_line_style( int s ){
    m_canvas->set_line_style(s);}

inline void Projector::set_col_index(int index){
    m_canvas->set_col_index(index);}

inline void Projector::set_cap_style(int style){
    m_canvas->set_cap_style(style);}

inline void Projector::set_marker_size(float size){
    m_canvas->set_marker_size(size);}

inline void Projector::set_enhanced(int on){
    m_canvas->set_enhanced(on);}

} // namespace gui

#endif //PROJECTOR_H

