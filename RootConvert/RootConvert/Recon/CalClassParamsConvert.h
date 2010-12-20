#ifndef RootConvert_CalClassParamsConvert_H
#define RootConvert_CalClassParamsConvert_H 1

#include <Event/Recon/CalRecon/CalClassParams.h>   
#include <reconRootData/CalClassParams.h>

namespace RootPersistence {

void convert( const Event::CalClassParams &, CalClassParams & ) ;
void convert( const CalClassParams &, Event::CalClassParams & ) ;

}

#endif
