//     $Id$
//  Author: Hardy, 29 Mar 94, T. Burnett after

#ifndef __PSDRAW_H
#define __PSDRAW_H


#include <fstream>
#include "gui/Draw2D.h"

namespace gui {

class PSdraw : public Draw2D
{
    
public:
    PSdraw (const char *filename="temp.ps", const char *format="DinA4");
    ~PSdraw();
    bool operator()(){return psfile.good();}
    // allow testing that file is open
    
    void draw_string (float x, float y, const char * title, int size=0);
    // draw a string at giveconst char *, int size= 0);
    
    void draw_string (const char* string, int size);
    // draw a string at current point
    
    void move_to (float, float);
    // move current point
    
    void line_to (float, float);
    // draw line from current point
    
    void draw_marker (float,float);
    // draw a marker at the point;
    
    void flush ();
    // flush drawing commands
    
    void set_quad (int quad, const char * title = 0, int selected=0);
    // define and clear viewport quadrant, optionally draw title and box
    
    void set_defaults();
    void set_line_style(int);
    struct PageFormat;
    
    virtual void set_rgb(float r,float g,float b);
    virtual void set_col_index(int);

private:
    std::ofstream  psfile;
    unsigned long  myforeground, mybackground;
    int            lastx, lasty; 
    int            screenWidth, screenHeight;
    bool           isStrokePending;
    int            flSingle, selected;
    
    struct HRectangle
    {
        int x, y, dx, dy;
        //HRectangle(int _x, int _y, int _dx, int _dy)   { x=_x; y=_y; dx=_dx; dy=_dy; };
    };
    HRectangle currentQuad;
    HRectangle screen;
    HRectangle quadrant[4];
    
    void resize (int xsize, int ysize);
    
    int xwin (float);
    int ywin (float);
    // private functions that convert to screen coordinates
    
        
};
    
} // namespace gui
#endif


