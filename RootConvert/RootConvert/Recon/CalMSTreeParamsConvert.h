#ifndef RootConvert_CalMSTreeParamsConvert_H
#define RootConvert_CalMSTreeParamsConvert_H 1

#include <Event/Recon/CalRecon/CalMSTreeParams.h>   
#include <reconRootData/CalMSTreeParams.h>

namespace RootPersistence {

void convert( const Event::CalMSTreeParams &, CalMSTreeParams & ) ;
void convert( const CalMSTreeParams &, Event::CalMSTreeParams & ) ;

}

#endif

