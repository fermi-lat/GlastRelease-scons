#ifndef _GlastApps_HelloWorld_H
#define _GlastApps_HelloWorld_H 1

#include "GaudiKernel/Algorithm.h"

class IGlastDetSvc;
class IGlastIRFLoadSvc;

/** @class CreateEvent.cpp
* @brief Test algorithm to insure that the TDS is filled appropriately with MCACDHit, 
MCCalorimeterHit and MCTKRHit data.

  Original Author:  Sawyer Gillespie

  To Run:
  - From VCMT
    -# Choose GlastSvc as the active package
    -# Choose test to be the active project
    -# Under "app parameters" type CreatEvent
    -# Click on the run button
  - From VC++
    -# Open the GlastSvc workspace
    -# Make test the active project and compile if necessary
    -# Update the project settings for test, setting the working directory under Debug to ".."
    -# Hit "F5" or click on the Build menu and choose Debug
* $Header$
*/

class CreateEvent : public Algorithm {

public:
  CreateEvent(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

  StatusCode testMcClass();

private:
    // the GlastDetSvc objects
    IGlastDetSvc*       m_detSvc;
    IGlastIRFLoadSvc*   m_irfLoadSvc;
};


#endif // _GlastApps_HelloWorld_H
