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
: CompositeMedium( (Medium*)0, new Box(size,size,size), "vacuum", (Detector*)0)
{
    setTitle("World");
    s_instance = this;
}

Medium&
World::addMedium(Medium* nextMedium)
{
    // intercept CompositeMedium to make sure that we are contain
    const Shape* pvol = &nextMedium->volume();
    if( pvol ) {

	float d = pvol->getMaxDimension();
	Point Q  = pvol->center();
	float newsize = sqrt(1.3333*(Q.mag2() + sqr(d)) );
	if( newsize > volume().getMaxDimension() ) {
	    delete _volume;
	    _volume = new Box(newsize,newsize,newsize);
	}
    }
    return CompositeMedium::addMedium(nextMedium);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
void World::createDetectorView(gui::DisplayRep& v)
{
    for(iterator it=begin(); it !=end(); ++it)
	(*it)->createDetectorView(v);
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


