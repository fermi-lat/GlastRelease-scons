// $Header$
//  Author: Toby Burnett
// Project: gui graphical display
//
// Defines a Display Representation that is tuned to the gui environment
// Almost all of the functionallity is derived from DisplayList. The most important are:
//
// 	void move_to(const GraphicsVector&);
//	// move current point to location
//
//	void line_to(const GraphicsVector&);
//	// draw a line to the specified point
//
//	void drawText(const std::string&);
//	// draw the text at the current point
//
//	void markerAt(const GraphicsVector&);
//	// draw a marker at the given point
//
// where GraphicsVector is a 3-vector with automatic conversion from SpacePoint
//
// This is an abstract class: an subclass *must* implement the method
//     virtual void update();
// which should create the display list by calls to line_to, etc. Another important virtual method is
//     virtual void clear();
// The default behaviour is to erase the display list, which may be bypassed by defining
// a dummy.
//

#ifndef DISPLAYREP_H
#define DISPLAYREP_H

#include "gui/DisplayList.h"
#include "gui/GraphicsVector.h"

namespace gui {

class DisplayRep :  public DisplayList
{
public:	
	// Implicit constructor
	// Implicit copy constructor
	// Implicit assignment operator
	// Implicit destructor
    enum LineStyle {	  /* Line styles for get/setlinestyle */
	SOLID_LINE= 0, DOTTED_LINE = 1,  CENTER_LINE= 2,  DASHED_LINE = 3, USERBIT_LINE = 4};


    void set_color(const std::string& color_name);
    // set the color for subsequent operations


    virtual void clear();
    // remove the representation

    void append(DisplayRep& other);
    void append(const DisplayRep& other);
    // append another representation to this one (removes representation)

    DisplayRep& nested();
    // return a representation that is nested within this one

    bool empty()const;
    // is it empty?

    struct ColorInfo { const char * name; unsigned r, g,b; };
    static ColorInfo pallete[];
    // color palette management
private:

    class Dummy;

};

} // namespace gui


#endif //DISPLAYREP_H

