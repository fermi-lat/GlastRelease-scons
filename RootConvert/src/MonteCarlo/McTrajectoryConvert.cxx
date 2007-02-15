#include <RootConvert/MonteCarlo/McTrajectoryConvert.h>
#include <RootConvert/Utilities/Toolkit.h>
#include "TRandom.h"

namespace RootPersistence {

void convert( const Event::McTrajectory & tdsTrajectory, McTrajectory & rootTrajectory ) 
{
    // Initialize data members
    rootTrajectory.initialize();
    rootTrajectory.Clear();

    return;
}
 
void convert( const McTrajectory & rootTrajectory, Event::McTrajectory & tdsTrajectory ) 
{
    return;
}

void convert( const Event::McTrajectoryPoint & tdsPoint, McTrajectoryPoint & rootPoint ) 
{
    rootPoint.initialize(convert(tdsPoint.getVolumeID()),
                          tdsPoint.getEnergy(),
                          convert3vector(tdsPoint.getPoint()));

    return;
}
 
void convert( const McTrajectoryPoint & rootPoint, Event::McTrajectoryPoint & tdsPoint ) 
{
    tdsPoint.setVolumeID(convert(rootPoint.getVolumeID()));
    tdsPoint.setEnergy(rootPoint.getEnergy());
    tdsPoint.setPoint(convert3vector(rootPoint.getPoint()));
        
    return;
}

}
