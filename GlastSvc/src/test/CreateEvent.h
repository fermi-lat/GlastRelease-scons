#ifndef _GlastApps_HelloWorld_H
#define _GlastApps_HelloWorld_H 1

#include "GaudiKernel/Algorithm.h"

class IGlastDetSvc;
class IGlastIRFLoadSvc;

/** @class CreateEvent.cpp
* @brief Test algorithm 
* 
* @author: Sawyer Gillespie
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
};


#endif // _GlastApps_HelloWorld_H
