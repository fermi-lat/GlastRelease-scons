#ifndef RootConvert_CalMomParamsConvert_H
#define RootConvert_CalMomParamsConvert_H 1

#include <Event/Recon/CalRecon/CalMomParams.h>   
#include <reconRootData/CalMomParams.h>

namespace RootPersistence {

void convert( const Event::CalMomParams &, CalMomParams & ) ;
void convert( const CalMomParams &, Event::CalMomParams & ) ;

}

#endif












