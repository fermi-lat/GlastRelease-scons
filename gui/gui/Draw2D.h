//     $Id$
//  Author: Toby Burnett
//
//
#ifndef DRAW2D_H
#define DRAW2D_H

namespace gui {

class Draw2D {

// Abstract base class for 2-dimensional event and detector plotting
// Functions as a canvas
public:
 
    virtual ~Draw2D(){};
    // destructor;
  

//--------------these methods must be defined by derived class-----

    virtual void draw_string(float x, float y, const char *, int=0)=0;
    // draw a string at given coords

    virtual void draw_string(const char*, int=0)=0;
    // draw a character string at current point

    virtual void move_to(float, float)=0;
    virtual void line_to(float, float)=0;
    // move or draw line from current point to specified one

    virtual void draw_marker(float,float)=0;
    // draw a marker (cross) at the point

    virtual void set_quad(int quad=0, const char * title = 0, int selected=0)=0;
    // define and clear viewport quadrant, optionally draw title and box
    // quad 0 means whole display area, 1-4 are respective quadrants of display
    // the drawer is expected to clip to the defined region

    enum Shading { NONE, BRIGHT, MEDIUM, DARK };
    virtual void fill_polygon(const float* xy, int npoints, Shading pattern);
    // draw and close the polygon, and fill it according to the shading. Implement below
    // as a simple draw

    //------------- following optionally implemented---------------------
    virtual void set_enhanced(int){};
    // turn on/off (device-dependent) enhancement attribute

    virtual void new_path(){}
    virtual void close_path(){}
    virtual void stroke() {}
    virtual void fill() {}
    // to set a newpath and stroke

    virtual void set_defaults(){};
    // set all attributes to default values

    virtual void flush(){};
    // flush the accumulated drawing commands


    virtual void set_line_width(float){};
    // set the line width

    virtual void set_col_index(int){};
    // set the color Index

    virtual void set_rgb(float,float,float){};
    virtual void set_line_style(int){};

    virtual void set_marker_size(float){};
    // change the current marker size

    virtual void set_cap_style(int){};


protected:

    Draw2D(){};
    // constructor is projected, since abstract class


private:

};
// define the default polygon, with no shading
inline void Draw2D::fill_polygon(const float* xy, int npoints, Draw2D::Shading )
{
    float x1 = *xy++, y1=*xy++;
    move_to(x1,y1);
    for(int i = 2; i<npoints; i+=2){
	float x = *xy++, y = *xy++; 
	line_to(x,y);
    }
    move_to(x1,y1);
}

} // namespace gui
#endif

