//     $Id$
//  Author: T. Burnett, G. Barrand
// Project: Event Display
//
// Generate 2-D drawing to an X11 window

#ifndef MOTIFGUI__XDRAW_H
#define MOTIFGUI__XDRAW_H

#include "gui/Draw2D.h"

#include <X11/Xlib.h>
using namespace gui;

class Xdraw : public Draw2D
{
 public:
   // ----------Interface to instantiator----------------

   Xdraw(int xsize=500, int ysize=500);
   // constructor
   ~Xdraw();

   void setDisplay(Display*,Window);
   // call to set display window, etc. Must be called before drawing

   class Rectangle
   {
     public: float x,y,dx,dy;
   }; // this is the same as something in X. XRectangle???

   void setXorMode(int);
   void drawRect(const Rectangle &);
   // set XOR-mode, draw selection rectangle

   void resize (int xsize=0, int ysize=0);
   int xwin(float)const;
   int ywin(float)const;
   bool clip(float x, float y, int& x1, int& y1, int& x2, int& y2); 
   // functions that convert to screen coordinates

   float xnorm(int)const;
   float ynorm(int)const;
   // functions that convert to normalized coordinates
   //

   // -----------User interface (via View3D)-------------
   void draw_string(float x, float y, const char * title, int size=0);
   // draw a string at giveconst char *, int size= 0);

   void draw_string(const char* string, int size);
   // draw a string at current point

   void move_to(float, float);
   // move current point

   void line_to(float, float);
   // draw line from current point

   void draw_marker(float,float);
   // draw a marker at the point;

   virtual void fill_polygon(const float* xy, int npoints, Shading pattern);

   void flush();

   void set_quad(int quad, const char * title = 0, int selected=0);
   // define and clear viewport quadrant, optionally draw title and box

   void set_defaults();
   void set_line_style(int);
   void set_col_index(int);

 private:
   Display       *mydisplay;
   Window         mywindow;
   int            myscreen;
   float          lastx, lasty;
   int            flSingle, selected;
   GC		gc;  // graphics context
   class HRectangle
      {
    public:
      int x, y, dx, dy;
      //HRectangle(int _x, int _y, int _dx, int _dy)   { x=_x; y=_y; dx=_dx; dy=_dy; };
      };

   HRectangle currentQuad;
   HRectangle screen;
   HRectangle quadrant[4];


};

#endif

