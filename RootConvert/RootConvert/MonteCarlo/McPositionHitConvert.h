
#ifndef RootConvert_McPositionHitConvert_H
#define RootConvert_McPositionHitConvert_H 1

#include <Event/MonteCarlo/McPositionHit.h>   
#include <mcRootData/McPositionHit.h>
#include "Rtypes.h"

namespace RootPersistence {

void convert( const Event::McPositionHit &, McPositionHit & ) ;
void convert( const McPositionHit &, Event::McPositionHit & ) ;

// for tests
void random( McPositionHit &, UInt_t i, Float_t randNum ) ;
bool compare( const McPositionHit &, const McPositionHit & ) ;

}

#endif












