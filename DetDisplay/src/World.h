// $Header$


#ifndef WORLD_H
#define WORLD_H

#include "CompositeMedium.h"

namespace gui {class DisplayRep; class Command;}
class Generator;

/// World: a special subclass of Medium intended to be the mother volume,
/// or "known world"
class World : public CompositeMedium
{
public:
    
    World(float size=1000);
    
    /// override from CompositeMedium to allow setting of box size
    virtual Medium&  addMedium(Medium* nextMedium);
        
    /// override in order to avoid displaying surrounding box
    void createDetectorView(gui::DisplayRep& v);
    /// return instance of special gui::DisplayRep objects
    gui::DisplayRep* detectorViewer();
    
    static World* instance();
    
private:
    static World* s_instance;
    class DetectorView;
};
#endif

