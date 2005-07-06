#ifndef RootConvert_CalClusterConvert_H
#define RootConvert_CalClusterConvert_H 1

#include <Event/Recon/CalRecon/CalCluster.h>   
#include <reconRootData/CalCluster.h>

namespace RootPersistence {

void convert( const Event::CalCluster &, CalCluster & ) ;
void convert( const CalCluster &, Event::CalCluster & ) ;

// for tests
void random( CalCluster &, UInt_t i, Float_t randNum ) ;
bool compare( const CalParams &, const CalParams & ) ;
bool compare( const CalClusterLayerData &, const CalClusterLayerData & ) ;
bool compare( const CalCluster &, const CalCluster & ) ;

}

#endif












