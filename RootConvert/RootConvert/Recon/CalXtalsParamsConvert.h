#ifndef RootConvert_CalXtalsParamsConvert_H
#define RootConvert_CalXtalsParamsConvert_H 1

#include <Event/Recon/CalRecon/CalXtalsParams.h>   
#include <reconRootData/CalXtalsParams.h>

namespace RootPersistence {

void convert( const Event::CalXtalsParams &, CalXtalsParams & ) ;
void convert( const CalXtalsParams &, Event::CalXtalsParams & ) ;

}

#endif
