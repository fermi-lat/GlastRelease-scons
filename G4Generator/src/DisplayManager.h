// $Header$

#ifndef DisplayManager_h
#define DisplayManager_h

// forward declarations
namespace gui {class DisplayControl; class DisplayRep; }

class HepTransform3D;
class Hep3Vector;
class Box;
#include <vector>


/**
    A simple class to manage a GuiSvc display of G4 objects. It is a singleton
  */
class DisplayManager {

public:
    //! ctor, with pointer to a DisplayControl object
    DisplayManager(gui::DisplayControl*);

    //! add a hit detector box to the display
    //! @param T global to local transformation
    //! @param x,y,z dimensions of the box
    void addBox(const HepTransform3D& T, double x, double y, double z);

    //! add to a the static display of all detector boxes
    void addDetectorBox(const Box* box);

    //! add a hit to the display 
    //! @param a,b initial, final points for the step
    void addHit(const Hep3Vector& a, const Hep3Vector& b);

    //! Add rep to display the tracks
    typedef std::vector<Hep3Vector> PointList;
    void addTrack(const PointList& track,  int charge);

    //! this is a singleton: access to the instance in this memory space
    static DisplayManager * instance(){return s_instance; }

private:
    gui::DisplayControl* m_display;

    gui::DisplayRep* m_boxes;
    gui::DisplayRep* m_hits;
    gui::DisplayRep* m_tracks;
    gui::DisplayRep* m_all_boxes;

    static DisplayManager* s_instance;
};


#endif