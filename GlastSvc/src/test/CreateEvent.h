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

/*! Short test algorithm to test the functionality of MCACDHit, MCCalorimeterHit, and MCTRKHit
*/

class CreateEvent : public Algorithm {

public:
  //! Constructor of this form must be provided
  CreateEvent(const std::string& name, ISvcLocator* pSvcLocator); 

  //! mandatory
  StatusCode initialize();
  //! mandatory
  StatusCode execute();
  //! mandatory
  StatusCode finalize();

private:
    // the GlastDetSvc objects
    IGlastDetSvc*       m_detSvc;
    IGlastIRFLoadSvc*   m_irfLoadSvc;
};


#endif // _GlastApps_HelloWorld_H
