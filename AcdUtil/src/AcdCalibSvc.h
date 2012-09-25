#ifndef AcdCalibSvc_H
#define AcdCalibSvc_H
// $Header$

// LOCAL 
#include "AcdCalibMgr.h"
#include "AcdCalibSvcBase.h"
#include "AcdUtil/IAcdCalibSvc.h"

// GLAST 
#include "CalibSvc/ICalibPathSvc.h"


// EXTLIB
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IIncidentListener.h"

class AcdCalibMgr;

// STD

/**
 * @class AcdCalibSvc
 *
 * \brief Instatiates AcdUtil::IAcdCalibSvc interface for calibrations used in reconstruction
 *
 * This service is used by the reconstruction code.  In particular, AcdRecon/AcdPha2MipTool.
 *
 * handles:
 * - data storage/destruction
 * - communication with Gleam lower level services
 * - checking of data validity period  
 * - extraction of acd-specific constants out of generic data objects
 * - creation/caching of local meta-data objects where needed
 *
 * Job Options:
 * - CalibDataSvc        ["CalibDataSvc"] Name of the gaudi svc used to access calib DB.
 * - DefaultFlavor       ["ideal"]        Flavor of calibration to use.
 * - FlavorPed           [""]             To override pedestal flavor.
 * - FlavorGain          [""]             To override gain (aka MIP peak) flavor.
 * - FlavorHighRange     [""]             To override high range calibration flavor.
 * - FlavorCoherentNoise [""]             To override coherent noise calibration flavor.
 *
 * \author  Eric Charles (from Zachary Fewtrell's CalCalib stuff)
 *
 **/

class AcdCalibSvc : public Service, public AcdCalibSvcBase, virtual public AcdUtil::IAcdCalibSvc {

public:

  AcdCalibSvc(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize () {return StatusCode::SUCCESS;}

  /// queryInterface - for implementing a Service this is necessary
  StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);

  /// return the service type
  const InterfaceID&  AcdCalibSvc::type () const {return IID_IAcdCalibSvc;}

  virtual MsgStream makeMessageLog() const {
    return MsgStream(msgSvc(),name());
  }
  
  virtual StatusCode retrieveObject(const std::string& path, DataObject*& pObject) {
    StatusCode sc = m_dataProviderSvc->retrieveObject(path, pObject);
	if (pObject) sc = m_dataProviderSvc->updateObject(pObject);
	return sc;
  }

  const std::string getCalibPath(const ICalibPathSvc::CalibItem item, const std::string& flavor="") const;

protected:

  virtual StatusCode getCalibMgr(AcdCalibData::CALTYPE type,				 
				 AcdCalibMgr*& calibMgr) {
    return getCalibrationMgr(type, calibMgr);
  }

  void addCalibration(AcdCalibMgr* calibMgr,const std::string& flavorName);


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
  // following classes all share many properties w/ AcdCalibSvc as they are
  // sort of 'employees' of AcdCalibSvc.  easiest way to do it is to make them
  // friends
  friend class AcdCalibMgr;
};

#endif // AcdCalibSvc_H
