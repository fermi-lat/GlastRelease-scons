// $Header$

#ifndef DisplayManager_h
#define DisplayManager_h

#include "idents/VolumeIdentifier.h"

// forward declarations
namespace gui {class DisplayControl; class DisplayRep;}

class HepTransform3D;
class Hep3Vector;
#include <vector>
#include <map>
#include <string>

/** 
 * @class DisplayManager
 *
 * @brief A singleton for GUI interaction
 * 
 * A simple class to manage a GuiSvc display of G4 objects. It is a singleton
 * 
 * @author T. Burnett
 */
class DisplayManager {

public:
    //! ctor, with pointer to a DisplayControl object
    DisplayManager(gui::DisplayControl*);

    //! add a hit detector box to the display
    //! @param T global to local transformation
    //! @param x,y,z dimensions of the box
    void addHitBox(const HepTransform3D& T, double x, double y, double z);

    //! add a integrating detector box to the display
    //! @param T global to local transformation
    //! @param x,y,z dimensions of the box
    void addIntegratingBox(const HepTransform3D& T, 
                           double x, double y, double z);

    //! add to the static display of all detector boxes
    //! @param detname The name of the logical volume
    //! @param T global to local transformation
    //! @param x,y,z dimensions of the box
    void addDetectorBox(std::string detname, const HepTransform3D& T, 
                        double x, double y, double z);

    //! add a hit to the display 
    //! @param a,b initial, final points for the step
    void addHit(const Hep3Vector& a, const Hep3Vector& b);

    //! Add rep to display the tracks
    //! @param track a list of dots to connect
    //! @param charge tell if neurtal
    //! @param primary flag if primary track, to distinguish from neutrals?
    typedef std::vector<Hep3Vector> PointList;
    void addTrack(const PointList& track,  int charge, bool primary=false);
    
    //! Add to rep that displays ID's
    void addIdDisplay(const HepTransform3D& T, idents::VolumeIdentifier id);

    //! this is a singleton: access to the instance in this memory space
    static DisplayManager * instance(){return s_instance; }

private:
    gui::DisplayControl* m_display;

    //! access all DisplayRep pointers by name in this map
    std::map<std::string, gui::DisplayRep*> m_detmap;

    static DisplayManager* s_instance;
};


#endif
