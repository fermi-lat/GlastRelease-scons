// $Header$
//  Author: Toby Burnett
// Project: gui graphical display
//

#ifndef DISPLAYREP_H
#define DISPLAYREP_H

#include "gui/DisplayList.h"
#include "gui/GraphicsVector.h"

namespace gui {

/** Defines a Display Representation that is tuned to the gui environment
 Almost all of the functionallity is derived from DisplayList. The most important are:

 -	void moveTo(const GraphicsVector&);
	// move current point to location

-	void lineTo(const GraphicsVector&);
	// draw a line to the specified point

-	void drawText(const std::string&);
	// draw the text at the current point

-	void markerAt(const GraphicsVector&);
	// draw a marker at the given point
<br>
 where GraphicsVector is a 3-vector with automatic conversion from Hep3Vector

 <br>This is an abstract class: an subclass *must* implement the method
     virtual void update();
 which should create the display list by calls to line_to, etc. Another important virtual method is
     virtual void clear();
 The default behaviour is to erase the display list, which may be bypassed by defining
 a dummy.
*/

class DisplayRep :  public DisplayList
{
public:	
	// Implicit constructor
	// Implicit copy constructor
	// Implicit assignment operator
	// Implicit destructor
    enum LineStyle {	  /* Line styles for get/setlinestyle */
	SOLID_LINE= 0, DOTTED_LINE = 1,  CENTER_LINE= 2,  DASHED_LINE = 3, USERBIT_LINE = 4};


    /// set the color for subsequent operations
    void setColor(const std::string& color_name);


    /// remove the representation
    virtual void clear();

    /// append another representation to this one (removes representation)
    void append(DisplayRep& other);
    void append(const DisplayRep& other);

    /// return a representation that is nested within this one
    DisplayRep& nested();

    /// is it empty?
    bool empty()const;

    /// color palette management
    struct ColorInfo { const char * name; unsigned r, g,b; };
    static ColorInfo pallete[];


    //! for compatability
    void line_to(const GraphicsVector& v){lineTo(v);}
    void move_to(const GraphicsVector& v){moveTo(v);}
    void set_color(const std::string& color_name){setColor(color_name);}
private:

    class Dummy;

};

} // namespace gui


#endif //DISPLAYREP_H

