#ifndef AcdCalibSvcBase_H
#define AcdCalibSvcBase_H
// $Header$

// LOCAL 
#include "AcdUtil/AcdCalib.h"

// GLAST 
#include "CalibSvc/ICalibPathSvc.h"

// EXTLIB
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/Property.h"
#include "GaudiKernel/MsgStream.h"

// STD
#include <map>
#include <vector>
#include <string>

// forward declares
class AcdCalibMgr;
class DataObject;

namespace CalibData {
  class AcdCalibObj;
}

/** @class AcdCalibSvcBase
    @author Eric Charles (from Zachary Fewtrell's CalCalib stuff)
    * \brief
    *
    * handles:
    * - data storage/destruction
    * - communication with Gleam lower level services
    * - checking of data validity period  
    * - extraction of acd-specific constants out of generic data objects
    * - creation/caching of local meta-data objects where needed
    *
    * \author  Eric Charles (from Zachary Fewtrell's CalCalib stuff)
    *
    */

class AcdCalibSvcBase : virtual public IIncidentListener {

public:

  AcdCalibSvcBase(); 
  virtual ~AcdCalibSvcBase();

  virtual const std::string getCalibPath(const ICalibPathSvc::CalibItem item, const std::string& flavor="") const = 0;

  virtual MsgStream makeMessageLog() const = 0;

  virtual StatusCode retrieveObject(const std::string& path, DataObject*& pObject) = 0;

  StatusCode getCalibrationMgr(AcdCalibData::CALTYPE type, AcdCalibMgr*& calib);


protected:

  /// Add a manager (and the associated flavor)
  void addMgr(AcdCalibData::CALTYPE type, AcdCalibMgr* mgr, const StringProperty* flavor);

  /// Add this to the Incident SVC
  StatusCode prepapreManagers(MsgStream& log, const std::string& defaultFlavor);

  /// hook the BeginEvent so that we can check our validity once per event.
  void handle ( const Incident& inc );

private:
  ////////////////////////////////////////////////
  ////// PARAMETER MANAGEMENT ////////////////////
  ////////////////////////////////////////////////


  //-- FRIEND CLASSES --//
  // following classes all share many properties w/ AcdSimCalibSvc as they are
  // sort of 'employees' of AcdSimCalibSvc.  easiest way to do it is to make them
  // friends
  friend class AcdCalibMgr;

  std::vector<std::pair<AcdCalibMgr*, const StringProperty*> >  m_mgrs;

};

#endif // AcdSimCalibSvc_H
