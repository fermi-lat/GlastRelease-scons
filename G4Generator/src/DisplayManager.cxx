// File and Version Information:
// $Header$
//
// Description: this class is responsable of the drawing of hits, hits volumes
// and trajectories if the GUI has been activated
//      
// Author(s):
//      T.Burnett

#include "DisplayManager.h"

#include "gui/DisplayControl.h"
#include "gui/SubMenu.h"

#include "geometry/Box.h"
#include "geometry/CoordTransform.h"
#include "geomrep/BoxRep.h"

#include "CLHEP/Geometry/Transform3D.h"
#include <map>
#include <cassert>

// this is the singleton static pointer
DisplayManager* DisplayManager::s_instance=0;

DisplayManager::DisplayManager( gui::DisplayControl* d)
:m_display(d)
{
    s_instance = this;

    //! minimal rep just to append stuff to. Allows default color
    class EmptyRep : public gui::DisplayRep {
    public:
        EmptyRep(std::string color="black"):m_color(color){}
        void update(){}
        void clear(){DisplayRep::clear(); setColor(m_color);}
    private: 
        std::string m_color;
    };

    d->add(m_detmap["steps"] = new EmptyRep, "hits", false);
    
    d->add(m_detmap["hit_boxes"] = new EmptyRep("orange"), 
           "hit Pos detectors");

    d->add(m_detmap["integrating_boxes"] = new EmptyRep("blue"), 
           "hit Int detectors");

    d->add((m_detmap["tracks"]= new EmptyRep("black")), "tracks");
    d->add((m_detmap["ids"] = new EmptyRep("black")), 
           "volume identifiers", false);
}
void DisplayManager::addDetectorBox(std::string detName, 
                                    const HepTransform3D& T, 
                                    double x, double y, double z)
{
    Box b(x,y,z);
    b.transform(CoordTransform(T.getRotation(), T.getTranslation()));
    gui::DisplayRep* rep = m_detmap[detName.substr(0,3)];
    assert(rep);
    rep->append(BoxRep(b));
}

void DisplayManager::addHitBox(const HepTransform3D& T,
                               double x, double y, double z)
{
    Box b(x,y,z);
    b.transform(CoordTransform(T.getRotation(), T.getTranslation()));
    //b.transform(T); 
    m_detmap["hit_boxes"]->append(BoxRep(b));
    
}
void DisplayManager::addIntegratingBox(const HepTransform3D& T,
                               double x, double y, double z)
{
    Box b(x,y,z);
    b.transform(CoordTransform(T.getRotation(), T.getTranslation()));
    m_detmap["integrating_boxes"]->append(BoxRep(b));
    
}
void DisplayManager::addHit( const Hep3Vector& a, const Hep3Vector& b)
{
    class LineRep : public gui::DisplayRep {
    public:
        LineRep(const Hep3Vector& a, const Hep3Vector& b) 
        {
            markerAt(a); moveTo(a);  lineTo(b);
        }

        void update(){}
    };

    m_detmap["steps"]->append( LineRep(a,b) );
}

void DisplayManager::addTrack(const PointList & track, int charge)
{
    class TrackRep : public gui::DisplayRep {
    public:
        TrackRep( const DisplayManager::PointList& track, int charge){
            setColor(charge==0? "white" : "black");
            DisplayManager::PointList::const_iterator pit = track.begin();
        moveTo(*pit++);
        for(; pit !=track.end(); ++pit) lineTo(*pit);
        }
        void update(){}
    };

    m_detmap["tracks"]->append(TrackRep(track,charge));
}
void DisplayManager::addIdDisplay(const HepTransform3D& T, 
                                  idents::VolumeIdentifier id)
{
    class IdRep : public gui::DisplayRep {
    public:
        IdRep(HepPoint3D center, std::string idname){
            moveTo( center);
            drawText(idname);
        }
        void update(){};
    };

    m_detmap["ids"]->append(IdRep(T.getTranslation(), id.name()) );
}
