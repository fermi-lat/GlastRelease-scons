// MaterialDisplay.h: interface for the MaterialDisplay class.
//
// $Header$
//  Author: T. Burnett
//////////////////////////////////////////////////////////////////////

#if !defined(MATERIALDISPLAY__INCLUDED_)
#define MATERIALDISPLAY__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "gui/DisplayRep.h"

namespace gui { class GuiMgr; class DisplayControl; class SubMenu;}
class Medium;
class Material;
class DisplayGeometry;

/** This is a special class that traverses the Medium hierarchy to generate displays
 of volumes according to the material.

 The constructor adds a "Materials" submenu to the Display menu. 
 This new submenu has a button to set the color for the next material display that is enabled.
 There is a button labelled with the name of each known material, which are initially hidden. 
 Clicking one generates the display, which creates the representations of all the volumes 
 that are used. 

 Note that if a medium is composite, the material only applies to the 
 volume that is not occupied by inner media. (So it should perhaps display also the next level?)
*/  
class MaterialDisplay 
{
public:
    MaterialDisplay(gui::GuiMgr*, DisplayGeometry* geom, gui::SubMenu& detmenu ); // ctor sets it all up
    ~MaterialDisplay(); // dtor to clean

    // nested representation class
    class Rep;
private:
    gui::DisplayControl* m_display;
    class ChooseColor; // GUI class
};

#endif // MATERIALDISPLAY__INCLUDED_


