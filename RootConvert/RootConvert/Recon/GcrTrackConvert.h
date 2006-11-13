#ifndef RootConvert_GcrTrackConvert_H
#define RootConvert_GcrTrackConvert_H 1


#include <Event/Recon/CalRecon/GcrReconClasses.h>   
#include <reconRootData/GcrTrack.h>

namespace RootPersistence {

void convert( const Event::GcrTrack &, GcrTrack & ) ;
void convert( const GcrTrack &, Event::GcrTrack & ) ;

}

#endif












