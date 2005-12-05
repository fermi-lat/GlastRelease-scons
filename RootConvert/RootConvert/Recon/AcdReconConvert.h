#ifndef RootConvert_AcdReconConvert_H
#define RootConvert_AcdReconConvert_H 1

#include <Event/Recon/AcdRecon/AcdRecon.h>   
#include <reconRootData/AcdRecon.h>

namespace RootPersistence {

void convert( const Event::AcdRecon &, AcdRecon & ) ;
void convert( const AcdRecon &, Event::AcdRecon & ) ;

}

#endif












