#include <RootConvert/MonteCarlo/McTkrStripConvert.h>
#include <RootConvert/Utilities/Toolkit.h>
#include "TRandom.h"

namespace RootPersistence {

void convert( const Event::McTkrStrip & tdsStrip, McTkrStrip & rootStrip ) 
{
    // Initialize data members
    rootStrip.Clear();

    const idents::VolumeIdentifier idTds = tdsStrip.getId();
    VolumeIdentifier idRoot;
    idRoot = RootPersistence::convert(idTds);

    rootStrip.initialize(idRoot, tdsStrip.getStripNumber(), tdsStrip.getEnergy(),
        tdsStrip.getNoise(), tdsStrip.getDeltaX(), tdsStrip.getDeltaY(), 0);

    return;
}
 
void convert( const McTkrStrip & rootstrip, Event::McTkrStrip & tdsstrip ) 
{
    return;
}

}
