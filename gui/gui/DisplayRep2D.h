//     $Id$
//  Author: T. Burnett
// Project: gui

#ifndef DisplayRep2D_H
#define DisplayRep2D_H 1


#include "gui/DisplayRep.h"
#include "gui/ViewPort.h"


namespace gui {

//Special Rep that draws directly to the canvas: subclass must implement
// the draw(Draw2D*) method
class DisplayRep2D : public DisplayRep {
public:
    
    virtual void draw(ViewPort* port);
    // override to draw directly to the canvas, accessible from the ViewPort
    
    virtual void draw2D(Draw2D* canvas)=0;
    // subclass must define to draw directly to the canvas
    
    virtual void update(){};
    // must define to be an DisplayRep: however, drawing is done with the draw above
    
private:
    
};

inline void DisplayRep2D::draw(ViewPort* port)
{
    if( !hidden() ) draw2D(port->canvas());
}

} // namespace gui

#endif //DisplayRep2D_H
