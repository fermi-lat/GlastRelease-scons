// Gui.cxx
// $Header$
//
// stuff to access the instance
#include "gui/GUI.h"

namespace gui {

GUI* GUI::s_instance=0;

GUI* GUI::instance(GUI* i){ 
    if(i) s_instance=i;
    return s_instance;
}

}