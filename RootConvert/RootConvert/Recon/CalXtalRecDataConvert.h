#ifndef RootConvert_CalXtalRecDataConvert_H
#define RootConvert_CalXtalRecDataConvert_H 1

#include <Event/Recon/CalRecon/CalXtalRecData.h>   
#include <reconRootData/CalXtalRecData.h>

namespace RootPersistence {

void convert( const Event::CalXtalRecData &, CalXtalRecData & ) ;
void convert( const CalXtalRecData &, Event::CalXtalRecData & ) ;

}

#endif












