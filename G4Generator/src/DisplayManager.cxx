// $Header$

#include "DisplayManager.h"

#include "gui/DisplayControl.h"

#include "geometry/Box.h"
#include "geometry/CoordTransform.h"
#include "geomrep/BoxRep.h"

#include "CLHEP/Geometry/Transform3D.h"



DisplayManager* DisplayManager::s_instance=0;

DisplayManager::DisplayManager( gui::DisplayControl* d)
:m_display(d)
{
    s_instance = this;
    
    class Boxes : public gui::DisplayRep {
    public:
        Boxes(){}
        void update(){ 
            setColor("blue");
        }
        void clear(){DisplayRep::clear(); setColor("blue");}
    };
    class AllBoxes : public gui::DisplayRep {
    public:
        AllBoxes(){setColor("grey");}
        void update(){}
        void clear(){}
    };
    class HitsRep : public gui::DisplayRep {
    public:
        HitsRep(){}
        void update(){
        }
        void clear(){DisplayRep::clear(); setColor("red");}
    private:
    };
    
    class TracksRep : public gui::DisplayRep {
    public:
        TracksRep(){}
        void update(){
        }
        void clear(){DisplayRep::clear(); setColor("black");
        }
    };
    d->add(m_all_boxes = new AllBoxes, "all detectors");

    d->add(m_hits= new HitsRep, "hits");
    
    d->add(m_boxes = new Boxes, "hit detectors");

    d->add(m_tracks = new TracksRep, "tracks");
}

void DisplayManager::addDetectorBox(const Box* b)
{
    m_all_boxes->append(BoxRep(*b));
}

void DisplayManager::addBox(const HepTransform3D& T, double x, double y, double z)
{
    Box b(x,y,z);
    b.transform(CoordTransform(T.getRotation(), T.getTranslation()));
    m_boxes->append(BoxRep(b));
    
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

    m_hits->append( LineRep(a,b) );
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

    m_tracks->append(TrackRep(track,charge));
}
