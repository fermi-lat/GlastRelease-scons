#ifndef RootConvert_CalFitParamsConvert_H
#define RootConvert_CalFitParamsConvert_H 1

#include <Event/Recon/CalRecon/CalFitParams.h>   
#include <reconRootData/CalFitParams.h>

namespace RootPersistence {

void convert( const Event::CalFitParams &, CalFitParams & ) ;
void convert( const CalFitParams &, Event::CalFitParams & ) ;

}

#endif












