#ifndef AcdCalibSvc_H
#define AcdCalibSvc_H
// $Header$

// LOCAL 
#include "AcdPedCalibMgr.h"
#include "AcdGainCalibMgr.h"

// GLAST 
#include "AcdUtil/IAcdCalibSvc.h"

// EXTLIB
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IIncidentListener.h"

// STD

/** @class AcdCalibSvc
    @author Eric Charles (from Zachary Fewtrell's CalCalib stuff)
    * \brief Instatiates IAcdCalibSvc interface, gets data from CalibDataSvc
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

class AcdCalibSvc : public Service, virtual public AcdUtil::IAcdCalibSvc, 
    virtual public IIncidentListener {

public:

  AcdCalibSvc(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize () {return StatusCode::SUCCESS;}

  /// queryInterface - for implementing a Service this is necessary
  StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);

  /// return the service type
  const InterfaceID&  AcdCalibSvc::type () const {return IID_IAcdCalibSvc;}

  /// get pedestal vals for given channel
  StatusCode getPedestal(idents::AcdId id, unsigned pmt,
			 CalibData::AcdPed*& pedestal){
    return m_pedMgr.getPed(id,pmt,pedestal);
  }

  /// get the mip peak vals for a given channel
  StatusCode getMipPeak(idents::AcdId id, unsigned pmt,
			CalibData::AcdGain*& mipPeak) {
    return m_gainMgr.getMipPeak(id,pmt,mipPeak);
  }

private:
  ////////////////////////////////////////////////
  ////// PARAMETER MANAGEMENT ////////////////////
  ////////////////////////////////////////////////

  // JobOptions PROPERTIES
  /// name of CalibDataSvc, main data source
  StringProperty m_calibDataSvcName;     

  ///  default flavor for all calib types, unless otherwise specified.
  StringProperty m_defaultFlavor;        
        
  /// calib flavor override for ped constants
  StringProperty m_flavorPed;            
  /// calib flavor override for mip peak constants
  StringProperty m_flavorMipPeak;            

  // GAUDI RESOURCES
  /// pointer to CalibDataSvc
  IService         *m_calibDataSvc;     

  /// pointer to IDataProviderSvc interface of CalibDataSvc
  IDataProviderSvc *m_dataProviderSvc;  
 

  /// The pedestal manager
  AcdPedCalibMgr    m_pedMgr;

  /// The gain (aka MIP peak) manager
  AcdGainCalibMgr   m_gainMgr;
  
  /// hook the BeginEvent so that we can check our validity once per event.
  void handle ( const Incident& inc );

  //-- FRIEND CLASSES --//
  // following classes all share many properties w/ AcdCalibSvc as they are
  // sort of 'employees' of AcdCalibSvc.  easiest way to do it is to make them
  // friends
  friend class AcdCalibMgr;
};

#endif // AcdCalibSvc_H
