#ifndef RootConvert_GcrXtalConvert_H
#define RootConvert_GcrXtalConvert_H 1

#include <Event/Recon/CalRecon/GcrReconClasses.h>   
#include <reconRootData/GcrXtal.h>

namespace RootPersistence {

void convert( const Event::GcrXtal &, GcrXtal & ) ;
void convert( const GcrXtal &, Event::GcrXtal & ) ;

}

#endif












