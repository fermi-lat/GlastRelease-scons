// CompositeDiffuse.h: interface for the CompositeDiffuse class.
//
//////////////////////////////////////////////////////////////////////
#include "FluxSvc/EventSource.h"
#include "CompositeSource.h"
//#include <vector>
#include <vector>

class FluxSource;

//! holds multiple Eventsource objects ; acts as a container for them.
class CompositeDiffuse: public CompositeSource { 
public:
    
    CompositeDiffuse(double totalFlux):
      m_totalFlux(totalFlux),m_unclaimedFlux(totalFlux)
      {};
      
      ~CompositeDiffuse(){}
      
      ///    add a source to the list
      //void addSource (EventSource* aSource);
      //void rmvSource (EventSource* aSource);
      
      /// generate an event from from one of the sources 
      /// which make up the composite, and return a pointer to it
      FluxSource* event (double time);
      
      /// Same as from CompositeSource, but allowing for total flux adjustment.
      void CompositeDiffuse::addSource (EventSource* aSource);
      
      ///Randomly determines a new source, and adds it.
      void addNewSource();
      
      double remainingFluxInterval();
      
private:
    double m_totalFlux; // The total flux from the entire sky
    double m_unclaimedFlux; // The amount of the flux "unaccounted for" by the known sources
    double getRandomFlux();
    long double pofi(long double intensity);
    long double logNlogS(long double flux);
    
    
};
