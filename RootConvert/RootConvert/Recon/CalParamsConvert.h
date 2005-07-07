#ifndef RootConvert_CalParamsConvert_H
#define RootConvert_CalParamsConvert_H 1

#include <Event/Recon/CalRecon/CalParams.h>   
#include <reconRootData/CalParams.h>

namespace RootPersistence {

void convert( const Event::CalParams &, CalParams & ) ;
void convert( const CalParams &, Event::CalParams & ) ;

}

#endif












