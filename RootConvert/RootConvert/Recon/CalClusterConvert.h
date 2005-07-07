#ifndef RootConvert_CalClusterConvert_H
#define RootConvert_CalClusterConvert_H 1

#include <Event/Recon/CalRecon/CalCluster.h>   
#include <reconRootData/CalCluster.h>

namespace RootPersistence {

void convert( const Event::CalCluster &, CalCluster & ) ;
void convert( const CalCluster &, Event::CalCluster & ) ;

}

#endif












