#ifndef AcdCalibSvc_H
#define AcdCalibSvc_H
// $Header$

// LOCAL 
#include "AcdPedCalibMgr.h"
#include "AcdGainCalibMgr.h"
#include "AcdVetoCalibMgr.h"
#include "AcdCnoCalibMgr.h"
#include "AcdRangeCalibMgr.h"
#include "AcdHighRangeCalibMgr.h"
#include "AcdCoherentNoiseCalibMgr.h"

// GLAST 
#include "AcdUtil/IAcdCalibSvc.h"
#include "CalibSvc/ICalibPathSvc.h"


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

  const std::string getCalibPath(const ICalibPathSvc::CalibItem item, const std::string& flavor="") const;

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

  /// get the veto vals for a given channel
  StatusCode getVeto(idents::AcdId id, unsigned pmt,
			CalibData::AcdVeto*& veto) {
    return m_vetoMgr.getVeto(id,pmt,veto);
  }

  /// get the cno vals for a given channel
  StatusCode getCno(idents::AcdId id, unsigned pmt,
			CalibData::AcdCno*& cno) {
    return m_cnoMgr.getCno(id,pmt,cno);
  }

  /// get the range crossover vals for a given channel
  StatusCode getRange(idents::AcdId id, unsigned pmt,
			CalibData::AcdRange*& range) {
    return m_rangeMgr.getRange(id,pmt,range);
  }

  /// get the high range vals for a given channel
  StatusCode geHighRange(idents::AcdId id, unsigned pmt,
			CalibData::AcdHighRange*& highRange) {
    return m_highRangeMgr.getHighRange(id,pmt,highRange);
  }

  /// get the high range vals for a given channel
  StatusCode getCoherentNoise(idents::AcdId id, unsigned pmt,
			      CalibData::AcdCoherentNoise*& calib) {
    return m_coherentNoiseMgr.getCoherentNoise(id,pmt,calib);
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

  ICalibPathSvc    *m_calibPathSvc;

  /// pointer to IDataProviderSvc interface of CalibDataSvc
  IDataProviderSvc *m_dataProviderSvc;   

  /// The pedestal manager
  AcdPedCalibMgr    m_pedMgr;

  /// The gain (aka MIP peak) manager
  AcdGainCalibMgr   m_gainMgr;

  /// The veto manager
  AcdVetoCalibMgr   m_vetoMgr;

  /// The cno manager
  AcdCnoCalibMgr   m_cnoMgr;

  /// The range (aka xover from low-high range) manager
  AcdRangeCalibMgr   m_rangeMgr;

  /// The high range manager
  AcdHighRangeCalibMgr   m_highRangeMgr;

  /// The coherent noise manager
  AcdCoherentNoiseCalibMgr   m_coherentNoiseMgr;


  /// hook the BeginEvent so that we can check our validity once per event.
  void handle ( const Incident& inc );

  //-- FRIEND CLASSES --//
  // following classes all share many properties w/ AcdCalibSvc as they are
  // sort of 'employees' of AcdCalibSvc.  easiest way to do it is to make them
  // friends
  friend class AcdCalibMgr;
};

#endif // AcdCalibSvc_H
