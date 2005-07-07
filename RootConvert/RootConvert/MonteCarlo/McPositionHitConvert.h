
#ifndef RootConvert_McPositionHitConvert_H
#define RootConvert_McPositionHitConvert_H 1

#include <Event/MonteCarlo/McPositionHit.h>   
#include <mcRootData/McPositionHit.h>
#include "Rtypes.h"

namespace RootPersistence {

void convert( const Event::McPositionHit &, McPositionHit & ) ;
void convert( const McPositionHit &, Event::McPositionHit & ) ;

}

#endif












