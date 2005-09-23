#ifndef RootConvert_CalMipTrackConvert_H
#define RootConvert_CalMipTrackConvert_H 1

#include <Event/Recon/CalRecon/CalMipClasses.h>   
#include <reconRootData/CalMipTrack.h>

namespace RootPersistence {

void convert( const Event::CalMipTrack &, CalMipTrack & ) ;
void convert( const CalMipTrack &, Event::CalMipTrack & ) ;

}

#endif












