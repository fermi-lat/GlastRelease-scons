#ifndef RootConvert_CalClusterConvert_H
#define RootConvert_CalClusterConvert_H 1

namespace Event { class CalCluster ; } // the TDS class
class CalCluster ;                     // the ROOT class

namespace RootPersistence {

void convert( const Event::CalCluster &, CalCluster & ) ;
void convert( const CalCluster &, Event::CalCluster & ) ;

}

#endif












