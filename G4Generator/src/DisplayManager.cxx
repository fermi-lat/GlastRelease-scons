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
    

    class AllBoxes : public gui::DisplayRep {
    public:
        AllBoxes(){setColor("grey");}
        void update(){}
        void clear(){}
    };

    class Boxes : public gui::DisplayRep {
    public:
        Boxes(){}
        void update(){}
        void clear(){DisplayRep::clear(); setColor("blue");}
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
    // key detector type with first 3 characters of name
    d->add(m_detmap["sid"]=m_detmap["top"]= new AllBoxes, "ACD");
    d->add(m_detmap["SiL"]= new AllBoxes, "TKR");
    d->add(m_detmap["CsI"]= new AllBoxes, "CAL");
    d->add(m_detmap["dio"]=new AllBoxes, "diodes", false);
    d->menu().addSeparator();

    d->add(m_detmap["steps"] = new HitsRep, "hits", false);
    
    d->add(m_detmap["hit_boxes"] = new Boxes, "hit detectors");

    d->add(m_detmap["tracks"]= new TracksRep, "tracks");
}
void DisplayManager::addDetectorBox(std::string detName, 
                                    const HepTransform3D& T, 
                                    double x, double y, double z)
{
    Box b(0.1*x,0.1*y,0.1*z);
    b.transform(CoordTransform(T.getRotation(), 0.1*T.getTranslation()));
    // b.transform(T);
    gui::DisplayRep* rep = m_detmap[detName.substr(0,3)];
    assert(rep);
    rep->append(BoxRep(b));
}

void DisplayManager::addHitBox(const HepTransform3D& T,
                               double x, double y, double z)
{
    Box b(0.1*x,0.1*y,0.1*z);
    b.transform(CoordTransform(T.getRotation(), 0.1*T.getTranslation()));
    //b.transform(T); 
    m_detmap["hit_boxes"]->append(BoxRep(b));
    
}
void DisplayManager::addHit( const Hep3Vector& a, const Hep3Vector& b)
{
    class LineRep : public gui::DisplayRep {
    public:
        LineRep(const Hep3Vector& a, const Hep3Vector& b) 
        {
            markerAt(0.1*a); moveTo(0.1*a);  lineTo(0.1*b);
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
        moveTo(0.1*(*pit++));
        for(; pit !=track.end(); ++pit) lineTo(0.1*(*pit));
        }
        void update(){}
    };

    m_detmap["tracks"]->append(TrackRep(track,charge));
}
