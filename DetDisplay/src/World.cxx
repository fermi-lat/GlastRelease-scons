// $Id$


#include "World.h"
#include "gui/DisplayRep.h"
#include "gui/Command.h"
#include "geometry/Box.h"
#include <cassert>
static inline double sqr(double x){return x*x;}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
World::World(float size)
: CompositeMedium( (Medium*)0, new Box(size,size,size), "vacuum")
{
    setTitle("World");
    s_instance = this;
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// implement singleton here
World* World::s_instance=0;

World* World::instance()
{
    if( s_instance==0 ) s_instance = new World;
    return s_instance;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Viewing support: derived classses, functions that return instance
// DetectorView assumes that detector never changes, doesn't allow deletion
class World::DetectorView : public gui::DisplayRep {
friend class World;
    void update();
    void clear(){} // avoid erasing
};

// can't not be inline because of statis variables
void World::DetectorView::update(){
    static bool set=false;
    if( set ) return;
    World::instance()->createDetectorView(*this);
    set = true;
}

gui::DisplayRep* World::detectorViewer()
{
    return new DetectorView;
}


