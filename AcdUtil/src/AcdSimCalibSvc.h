#ifndef AcdSimCalibSvc_H
#define AcdSimCalibSvc_H
// $Header$

// LOCAL 
#include "AcdCalibMgr.h"

// GLAST 
#include "AcdUtil/IAcdCalibSvc.h"
#include "CalibSvc/ICalibPathSvc.h"


// EXTLIB
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IIncidentListener.h"

// STD

/** @class AcdSimCalibSvc
    @author Eric Charles (from Zachary Fewtrell's CalCalib stuff)
    * \brief Instatiates IAcdSimCalibSvc interface, gets data from CalibDataSvc
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

class AcdSimCalibSvc : public Service, public AcdCalibSvcBase, virtual public AcdUtil::IAcdCalibSvc {

public:
  
  AcdSimCalibSvc(const std::string& name, ISvcLocator* pSvcLocator); 
  
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize () {return StatusCode::SUCCESS;}

  /// queryInterface - for implementing a Service this is necessary
  StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);

  /// return the service type
  const InterfaceID&  AcdSimCalibSvc::type () const {return IID_IAcdCalibSvc;}

  virtual MsgStream makeMessageLog() const {
    return MsgStream(msgSvc(),name());
  }
  
  virtual StatusCode retrieveObject(const std::string& path, DataObject*& pObject) {
    return m_dataProviderSvc->retrieveObject(path, pObject);
  }

  const std::string getCalibPath(const ICalibPathSvc::CalibItem item, const std::string& flavor="") const;

protected:

  virtual StatusCode getCalibMgr(AcdCalibData::CALTYPE type,				 
				 AcdCalibMgr*& calibMgr) {
    return getCalibrationMgr(type, calibMgr);
  }

  void addCalibration(AcdCalibMgr* calibMgr, const std::string& flavorName);

private:
  ////////////////////////////////////////////////
  ////// PARAMETER MANAGEMENT ////////////////////
  ////////////////////////////////////////////////

  // JobOptions PROPERTIES
  /// name of CalibDataSvc, main data source
  StringProperty m_calibDataSvcName;     

  ///  default flavor for all calib types, unless otherwise specified.
  StringProperty m_defaultFlavor;        
        

  // GAUDI RESOURCES
  /// pointer to CalibDataSvc
  IService         *m_calibDataSvc;   

  ICalibPathSvc    *m_calibPathSvc;

  /// pointer to IDataProviderSvc interface of CalibDataSvc
  IDataProviderSvc *m_dataProviderSvc;   

  //-- FRIEND CLASSES --//
  // following classes all share many properties w/ AcdSimCalibSvc as they are
  // sort of 'employees' of AcdSimCalibSvc.  easiest way to do it is to make them
  // friends
  friend class AcdCalibMgr;
};

#endif // AcdSimCalibSvc_H
