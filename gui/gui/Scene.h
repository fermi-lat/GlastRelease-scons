//     $Id$
//  Author: Toby Burnett
//

#ifndef SCENE_H
#define SCENE_H

#include "gui/DisplayRep.h"
#include <list>

namespace gui {
class ViewPort;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Scene 
{
public:
    Scene();
    ~Scene();
    
  
    void add(DisplayRep*);
    // add to the rep to the list to display 
	
    void remove(DisplayRep*);
    // remove from the list

    void draw(ViewPort *v);
    // have the members of the list draw themselves onto the ViewPort

    void clear();
    // send clear message to all attached reps

    void update();
    // ask all representations to update

    void enable();
    // enable the display

    void disable();
    // disable the display

    bool enabled()const;
    // return enabled status


private:
    typedef std::list<DisplayRep*> RepList;
    RepList m_replist;	// list of Reps
    
    bool    m_enabled;	// can disable display 

};
//////// inlines //////////////////
inline bool Scene::enabled()const{return m_enabled;}

} // namespaced gui

#endif //SCENE_H






