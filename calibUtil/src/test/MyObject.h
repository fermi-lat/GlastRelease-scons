// $Header$
/// file for sample client

#ifndef CALIBUTIL_MYOBJECT_H
#define CALIBUTIL_MYOBJECT_H

#include "calibUtil/ClientObject.h"

// Used by test_strips.cxx test program.  See that file for implementation.
class MyObject : public calibUtil::ClientObject {
public: 
    
  MyObject() {};
  virtual ~MyObject() {};
  /// Performs client specified function on the data, one list at a time
  calibUtil::StripSrv::eRet 
  readData(calibUtil::StripSrv::towerRC towerId, 
           unsigned int trayNum, 
           calibUtil::StripSrv::eUnilayer uni, 
           calibUtil::StripSrv::eBadness howBad,
           const calibUtil::StripSrv::StripCol* const strips);
};
#endif
  
