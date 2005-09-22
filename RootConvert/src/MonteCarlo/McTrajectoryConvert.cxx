#include <RootConvert/MonteCarlo/McTrajectoryConvert.h>
#include <RootConvert/Utilities/Toolkit.h>
#include "TRandom.h"

namespace RootPersistence {

void convert( const Event::McTrajectory & tdsTrajectory, McTrajectory & rootTrajectory ) 
{
    // not to be managed by RootConvert
    McParticle * mcPartRoot = 0 ;

    // Initialize data members
    rootTrajectory.initialize(mcPartRoot, tdsTrajectory.getCharge());

    // Retrieve points on trajectory
    const std::vector<Hep3Vector> tdsPoints = tdsTrajectory.getPoints();
    std::vector<Hep3Vector>::const_iterator tdsPointsIter;

    // Loop through hits creating root equivalents
    for(tdsPointsIter = tdsPoints.begin(); tdsPointsIter != tdsPoints.end(); tdsPointsIter++)
    {
        const Hep3Vector tdsPoint = *tdsPointsIter;
        TVector3* rootPoint = new TVector3(tdsPoint.x(),tdsPoint.y(),tdsPoint.z());

        rootTrajectory.addMcPoint(rootPoint);
    }

    return;
}
 
void convert( const McTrajectory & rootTrajectory, Event::McTrajectory & tdsTrajectory ) 
{
    // Not to be managed by RootConvert
    Event::McParticle* mcPartTds = 0;

    // Set this to make sure it is something
    tdsTrajectory.setMcParticle(mcPartTds);

    // Now for charge
    tdsTrajectory.setCharge(rootTrajectory.getMcTrajectoryCharge());

    // Retrieve points read back in
    TObjArray* rootPoints = rootTrajectory.getMcPointCol();

    // std vector to store points
    std::vector<Hep3Vector> tdsPoints;

    // Loop through root points and convert
    for(int idx = 0; idx < rootPoints->GetEntries(); idx++)
    {
        TVector3*  rootPoint = (TVector3 *)rootPoints->At(idx);
        
        tdsPoints.push_back(Hep3Vector(rootPoint->x(), rootPoint->y(), rootPoint->z()));
    }

    // Add these points to TDS object
    tdsTrajectory.addPoints(tdsPoints);

    return;
}

}
