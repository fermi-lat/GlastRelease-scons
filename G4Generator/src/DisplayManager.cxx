// $Header$

#include "DisplayManager.h"

#include "gui/DisplayControl.h"
#include "gui/SubMenu.h"

#include "geometry/Box.h"
#include "geometry/CoordTransform.h"
#include "geomrep/BoxRep.h"

#include "CLHEP/Geometry/Transform3D.h"
#include <map>
#include <cassert>


DisplayManager* DisplayManager::s_instance=0;

DisplayManager::DisplayManager( gui::DisplayControl* d)
:m_display(d)
{
    s_instance = this;
    

    class Boxes : public gui::DisplayRep {
    public:
        Boxes(std::string color="blue"):m_color(color){}
        void update(){}
        void clear(){DisplayRep::clear(); setColor(m_color);}
    private: 
        std::string m_color;
    };
    class HitsRep : public gui::DisplayRep {
    public:
        HitsRep(){}
        void update(){}
        void clear(){DisplayRep::clear(); setColor("red");}
    private:
    };
    
    class TracksRep : public gui::DisplayRep {
    public:
        TracksRep(){}
        void update(){}
        void clear(){DisplayRep::clear(); setColor("black"); }
    };
    d->add(m_detmap["steps"] = new HitsRep, "hits", false);
    
    d->add(m_detmap["hit_boxes"] = new Boxes("blue"), "hit detectors");

    d->add(m_detmap["integrating_boxes"] = new Boxes("orange"), "hit integrating detectors");

    d->add(m_detmap["tracks"]= new TracksRep, "tracks");
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
    //b.transform(T); 
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
