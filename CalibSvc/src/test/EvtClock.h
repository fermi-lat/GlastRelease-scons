//$Header$
#ifndef EVTCLOCK_H
#define EVTCLOCK_H 1

#include "facilities/Timestamp.h"

// Base class
#include "GaudiKernel/Algorithm.h"

// Forward declarations
class IDetDataSvc;

/** @class EvtClock EvtClock.h

    Simple algorithm to set fake event times on the detector data service.

    @author Joanne Bogart, adapted from similar class of Andrea Valassi 
*/

class EvtClock : public Algorithm {

 public:

  EvtClock( const std::string& name, ISvcLocator* pSvcLocator ); 
  
  // Algorithm standard methods
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();
  
 private:

  /// Absolute time of current event
  //  longlong i_evtTime();
  facilities::Timestamp i_evtTime();

 private:

  /// Current event number
  long m_eventNumber;

  /// Absolute time of first event (yyyy-mm-dd_hh:mm, trailing fields
  /// optional)
  std::string m_startTimeAsc;

  /// Absolute time of first event (seconds)
  long m_startTime;

  /// Absolute time spacing between events
  long m_delayTime;

  /// Handle to the IDetDataSvc interface of the CalibDataSvc
  IDetDataSvc* m_detDataSvc;
  
};

#endif    // EVTCLOCK_H
