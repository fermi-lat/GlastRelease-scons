
#ifndef RootConvert_McTrajectoryConvert_H
#define RootConvert_McTrajectoryConvert_H 1

#include <Event/MonteCarlo/McTrajectory.h>   
#include <mcRootData/McTrajectory.h>
#include "Rtypes.h"

namespace RootPersistence {

void convert( const Event::McTrajectory &, McTrajectory & ) ;
void convert( const McTrajectory &, Event::McTrajectory & ) ;

}

#endif












