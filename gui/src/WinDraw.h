//     $Header$
//
//
// Special Draw2D subclass for Win graphics
// -----------------------------
#ifndef __WinDRAW_H
#define __WinDRAW_H

#include "gui/Draw2D.h"
#define STRICT
#include <windows.h>

using namespace gui;

class WinDraw : public Draw2D {
    
public:
    WinDraw(HWND win);
    ~WinDraw();
    
    void draw_string(float x, float y, const char * title, int size=0);
    // draw a string at given coords
    
    void draw_string(const char *, int size= 0);
    // draw a string at current point
    
    void move_to(float, float);
    // move current point
    
    void line_to(float, float);
    // draw line from current point
    
    virtual void fill_polygon(const float* xy, int npoints, Shading pattern);
    // draw the polygon, fill it according to the shading

    void draw_marker(float,float);
    // draw a marker at the point;
    
    void set_quad(int quad, const char * title = 0, int selected=0);
    // define and clear viewport quadrant, optionally draw title and box
    
    void drawTitle(const char * title, int selected);
    // label current view
    
    void set_defaults();
    // set default drawing attributes
    
    void drawRectangle();
    // draw bounding rectangle
    
    void clearArea();
    // clear the display area
    
    void clearArea(int,int,int,int);
    // set viewport and clear it
    
    void flush(){};
    // flush the accumulated drawing commands
    
    
    void set_line_width(float);
    // set the line width
    
    void set_col_index(int in);
    // set the color Index
    
    void set_rgb(float,float,float);
    void set_line_style(int);
    
    void set_marker_size(float);
    
    void setHDC(HDC _hDC);
    // special to WinDraw: must be set before drawing

    float xnorm(int);
    float ynorm(int);
    // transform from window to 2-coords

    void invalidate();
    // force redraw by invalidating

private:
    int maxx, maxy;         // coords of viewport
    bool clip(float x, float y,
		 int& ix1, int& iy1, int& ix2, int& iy2);
    float m_lastx, m_lasty;
    HDC m_hdc;
    HWND m_win;
};
#endif

