// $Header$
 
#ifndef CALIBUTIL_CLIENTOBJECT_H
#define CALIBUTIL_CLIENTOBJECT_H

#include "calibUtil/StripSrv.h"

namespace calibUtil {

  class ClientObject {
  public: 
    
    /// Performs client specified function on data
    virtual unsigned int readData(towerRC towerId, unsigned int trayNum, 
                      uniL uniLayer, std::vector<unsigned int> v) = 0;

  };

}// end of namespace calibUtil

#endif


