// $Header$
#ifndef _GlastApps_HelloWorld_H
#define _GlastApps_HelloWorld_H 1


// Include files
#include "Gaudi/Algorithm/Algorithm.h"


//--------------------------------------------------------------------
//
//  CreateEvent: attempts to create a GlastEvent tree in Gaudi
//
//  Author :   Sawyer Gillespie
//
//--------------------------------------------------------------------


// forward declarations
class IGlastDetSvc;
class IGlastIRFLoadSvc;

/**
  Attempts to create a GlastEvent tree in Gaudi.

  @Author :   Sawyer Gillespie
 */

class CreateEvent : public Algorithm {

public:
  /// Constructor of this form must be provided
  CreateEvent(const std::string& name, ISvcLocator* pSvcLocator); 

  // Three mandatory member functions of any algorithm
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

private:
    // the GlastDetSvc objects
    IGlastDetSvc*       m_detSvc;
    IGlastIRFLoadSvc*   m_irfLoadSvc;
};


#endif // _GlastApps_HelloWorld_H
