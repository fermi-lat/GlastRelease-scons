// CompositeDiffuse.h: interface for the CompositeDiffuse class.
/** 
* \class CompositeDiffuse
*
* \brief CompositeDiffuse functions like CompositeSource, with added point sources
* CompositeDiffuse takes not only a list of initial sources, but also a total flux 
* over the whole sky.  It then attempts to fill in the remaining flux by generating point sources according
* to a logN/logS characteristic.
* \author Sean Robinson, University of Washington, 2002
* 
* $Header $
*/

#include "FluxSvc/EventSource.h"
#include "CompositeSource.h"
#include <vector>

class FluxSource;

//! holds multiple Eventsource objects ; acts as a container for them.
class CompositeDiffuse: public CompositeSource { 
public:
    
    CompositeDiffuse(double totalFlux):
      m_totalFlux(totalFlux),m_unclaimedFlux(totalFlux)
      {};
      
      ~CompositeDiffuse(){}
      
      
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
