//     $Id$
//  Author: Toby Burnett

#include "gui/DisplayRep.h"

namespace gui {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class DisplayRep::Dummy : public DisplayRep {
// Create and return a new "nested" DisplayRep object. This is only
// to allow the nested show/hide behaviour defined by DisplayList.
    friend class DisplayRep;
    void update(){}  
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DisplayRep&
DisplayRep::nested()
{
    DisplayRep& inner = *new Dummy;
    addDisplayList(&(inner));
    return inner;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void
DisplayRep::clear()
{
	deleteRepresentation();
}
bool
DisplayRep::empty()const
{
    return !hasRepresentation();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void
DisplayRep::append(DisplayRep& other)
{
    other.update();
    DisplayList::append(&other);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void
DisplayRep::append(const DisplayRep& other)
{
    append(const_cast<DisplayRep&>(other));
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void
DisplayRep::setColor(const std::string& color_name)
{
    int i=0;
    while(pallete[i].name ) {
        if( std::string(pallete[i].name) == color_name ){
            set_col_index(i);
            return;
        }
        i++; 
    }
    set_color("black"); // failed to find the name
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DisplayRep::ColorInfo 
DisplayRep::pallete[]= {
     { "white",    255, 255, 255},
     { "black",      0,   0,  0 },
     { "wheat",    245, 222, 179},
     { "violet",   238, 130, 238},
     { "turquoise", 64, 224, 208},
     { "thistle",  216, 191, 216},
     { "tan",      210, 180, 140},
     { "sienna",   160,  82,  45},
     { "salmon",   250, 128, 114},
     { "red",      255,   0,   0},
     { "plum",     221, 160, 221},
     { "pink",     255, 192, 203},
     { "orchid",   218, 112, 214},
     { "orange",   255, 165,   0},
     { "maroon",   176,  48,  96},
     { "magenta",  255,   0, 255},
     { "khaki",    240, 230, 140},
     { "grey",     192, 192, 192},
     { "green",      0, 255,   0},
     { "gold",     255, 215,   0},
     { "brown",    165,  42,  42},
     { "cyan",       0, 255, 255},
     { "coral",    255, 127,  80},
     { "navy",       0,   0, 128},
     { "blue",       0,   0, 255},
     { "yellow",   255, 255,   0},
     { "aquamarine",127, 255, 212},
     { 0, 0, 0,0} // tag end with null name
};

} // namespace gui
