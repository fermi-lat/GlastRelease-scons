
#ifndef RootConvert_McTkrStripConvert_H
#define RootConvert_McTkrStripConvert_H 1

#include <Event/MonteCarlo/McTkrStrip.h>   
#include <mcRootData/McTkrStrip.h>
#include "Rtypes.h"

namespace RootPersistence {

void convert( const Event::McTkrStrip &, McTkrStrip & ) ;
void convert( const McTkrStrip &, Event::McTkrStrip & ) ;

}

#endif












