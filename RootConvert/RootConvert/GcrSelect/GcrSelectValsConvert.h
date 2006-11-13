#ifndef RootConvert_GcrSelectValsConvert_H
#define RootConvert_GcrSelectValsConvert_H 1


#include <Event/Recon/CalRecon/GcrSelectClasses.h>   
#include <gcrSelectRootData/GcrSelectVals.h>

namespace RootPersistence {

void convert( const Event::GcrSelectVals &, GcrSelectVals & ) ;
void convert( const GcrSelectVals &, Event::GcrSelectVals & ) ;

}

#endif












