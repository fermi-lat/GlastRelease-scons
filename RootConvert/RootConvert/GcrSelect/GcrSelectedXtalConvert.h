#ifndef RootConvert_GcrSelectedXtalConvert_H
#define RootConvert_GcrSelectedXtalConvert_H 1

#include <Event/Recon/CalRecon/GcrSelectClasses.h>   
#include <gcrSelectRootData/GcrSelectedXtal.h>

namespace RootPersistence {

void convert( const Event::GcrSelectedXtal &, GcrSelectedXtal & ) ;
void convert( const GcrSelectedXtal &, Event::GcrSelectedXtal & ) ;

}

#endif












