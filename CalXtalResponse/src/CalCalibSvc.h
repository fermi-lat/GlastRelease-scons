#ifndef CalCalibSvc_H
#define CalCalibSvc_H

// LOCAL 
#include "AsymMgr.h"
#include "IntNonlinMgr.h"
#include "MPDMgr.h"
#include "PedMgr.h"
#include "TholdCIMgr.h"
#include "TholdMuonMgr.h"
#include "IdealCalCalib.h"

// GLAST 
#include "CalXtalResponse/ICalCalibSvc.h"

// EXTLIB
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IIncidentListener.h"

// STD

using namespace CalDefs;

/** @class CalCalibSvc
 * \brief Instatiates ICalCalibSvc interface, retrieves data from CalibDataSvc
 *
 * handles:
 * - data storage/destruction
 * - communication with Gleam lower level services
 * - checking of data validity period  
 * - extraction of cal-specific constants out of generic data objects
 * - creation/caching of spline function objects where needed.
 *
 * \author  Zachary Fewtrell
 *
 */

class CalCalibSvc : public Service, virtual public ICalCalibSvc, 
 virtual public IIncidentListener {

 public:

  CalCalibSvc(const string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize () {return StatusCode::SUCCESS;}



  /// queryInterface - for implementing a Service this is necessary
  StatusCode queryInterface(const IID& riid, void** ppvUnknown);

  static const InterfaceID& interfaceID() {return ICalCalibSvc::interfaceID();}

  /// return the service type
  const IID&  CalCalibSvc::type () const {return IID_ICalCalibSvc;}

  /// retrieve MeVPerDac ratios for given xtal
  StatusCode getMeVPerDac(const CalXtalId &xtalId,
                          CalibData::ValSig &asymLrg,
                          CalibData::ValSig &asymSm) {
    return m_mpdMgr.getMPD(xtalId, asymLrg, asymSm);
  }

  /// retrieve integral non-linearity vals for given xtal/face/rng
  StatusCode getIntNonlin(const CalXtalId &xtalId,
                          const vector< float > *&adcs,
                          const vector< unsigned > *&dacs,
                          float &error) {
    return m_intNonlinMgr.getIntNonlin(xtalId, adcs, dacs, error);
  }

  /// retrieve pedestal vals for given xtal/face/rng
  StatusCode getPed(const CalXtalId &xtalId,
                    float &avr,
                    float &sig,
                    float &cos) {
    return m_pedMgr.getPed(xtalId, avr, sig, cos);
  }

  /// retrieve Asymmetry calibration information for one xtal
  StatusCode getAsym(const CalXtalId &xtalId,
                     const vector<CalibData::ValSig> *&asymLrg,
                     const vector<CalibData::ValSig> *&asymSm,
                     const vector<CalibData::ValSig> *&AsymNSPB,
                     const vector<CalibData::ValSig> *&asymPSNB,
                     const vector<float> *&xVals) {
    return m_asymMgr.getAsym(xtalId, asymLrg, asymSm, AsymNSPB, asymPSNB, xVals);
  }

  /// retrieve threshold calibration constants as measured w/ charnge injection
  StatusCode getTholdCI(const CalXtalId &xtalId,
                        CalibData::ValSig &FLE,
                        CalibData::ValSig &FHE,
                        CalibData::ValSig &LAC
                        ) {
    return m_tholdCIMgr.getTholds(xtalId, FLE, FHE, LAC);
  }

  /// retrieve Upper Level Discriminator threshold as measured w/ charnge injection for given xtal/face/rng
  StatusCode getULDCI(const CalXtalId &xtalId,
                      CalibData::ValSig &ULDThold) {
    return m_tholdCIMgr.getULD(xtalId, ULDThold);
  }

  /// retrieve pedestal calibration constants as measured during charge injection threshold testing.
  StatusCode getPedCI(const CalXtalId &xtalId,
                      CalibData::ValSig &ped) {
    return m_tholdCIMgr.getPed(xtalId, ped);
  }

  /// retrieve threshold calibration constants as measured w/ muon calibration
  StatusCode getTholdMuon(const CalXtalId &xtalId,
                          CalibData::ValSig &FLE,
                          CalibData::ValSig &FHE
                          ) {
    return m_tholdMuonMgr.getTholds(xtalId, FLE, FHE);
  }

  /// retrieve pedestal calibration constants as measured during muon calibration threshold testing.
  StatusCode getPedMuon(const CalXtalId &xtalId,
                        CalibData::ValSig &ped) {
    return m_tholdMuonMgr.getPed(xtalId, ped);
  }

  StatusCode evalDAC(const CalXtalId &xtalId, double adc, double &dac) {
    return m_intNonlinMgr.evalDAC(xtalId, adc, dac);
  }
  
  StatusCode evalADC(const CalXtalId &xtalId, double dac, double &adc) {
    return m_intNonlinMgr.evalADC(xtalId, dac, adc);
  }
  
  StatusCode evalAsymLrg(const CalXtalId &xtalId, double pos, double &asymLrg) {
    return m_asymMgr.evalAsymLrg(xtalId, pos, asymLrg);
  }
  
  StatusCode evalPosLrg(const CalXtalId &xtalId, double asymLrg, double &pos) {
    return m_asymMgr.evalPosLrg(xtalId, asymLrg, pos);
  }
  
  StatusCode evalAsymSm(const CalXtalId &xtalId, double pos, double &asymSm) {
    return m_asymMgr.evalAsymSm(xtalId, pos, asymSm);
  }
  
  StatusCode evalPosSm(const CalXtalId &xtalId, double asymSm, double &pos) {
    return m_asymMgr.evalPosSm(xtalId, asymSm, pos);
  }
  
  StatusCode evalAsymNSPB(const CalXtalId &xtalId, double pos, double &asymNSPB) {
    return m_asymMgr.evalAsymNSPB(xtalId, pos, asymNSPB);
  }
  
  StatusCode evalPosNSPB(const CalXtalId &xtalId, double asymNSPB, double &pos) {
    return m_asymMgr.evalPosNSPB(xtalId, asymNSPB, pos);
  }
  
  StatusCode evalAsymPSNB(const CalXtalId &xtalId, double pos, double &asymPSNB) {
    return m_asymMgr.evalAsymPSNB(xtalId, pos, asymPSNB);
  }
  
  StatusCode evalPosPSNB(const CalXtalId &xtalId, double asymPSNB, double &pos) {
    return m_asymMgr.evalPosPSNB(xtalId, asymPSNB, pos);
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

  /// calib flavor override for int-nonlin constants
  StringProperty m_flavorIntNonlin;      
  /// calib flavor override for asymmetry constants
  StringProperty m_flavorAsym;           
  /// calib flavor override for ped constants
  StringProperty m_flavorPed;            
  /// calib flavor override for MeVPerDac constants
  StringProperty m_flavorMPD;            
  /// calib flavor override for CI measured thresholds
  StringProperty m_flavorTholdCI;        
  /// calib flavor override for Muon measure thresholds
  StringProperty m_flavorTholdMuon;      

  /// xml file contains 'ideal' flavor parameters
  StringProperty m_idealCalibXMLPath;

  /**
     \brief Skips time-consuming MsgStream object creations in frequently called functions
  
     Basically we cannot be creating millions MsgStream objects which are never used
     as would happen if the following code were in a frequently called function
  
     \code
     MsgStream msglog; 
     msglog << MSG::VERBOSE << "almost never used, but takes some time" << endreq;
     \endcode
  */
  BooleanProperty m_superVerbose;

  // GAUDI RESOURCES
  IService         *m_calibDataSvc;     ///< pointer to CalibDataSvc

  /// pointer to IDataProviderSvc interface of CalibDataSvc
  IDataProviderSvc *m_dataProviderSvc;  

  /// stores 'ideal' flavor generic calibration constants
  IdealCalCalib m_idealCalib;

  /// Load ideal calibration values from xml file
  StatusCode CalCalibSvc::loadIdealCalib() {
    return m_idealCalib.readCfgFile(m_idealCalibXMLPath);
  }

  PedMgr       m_pedMgr;
  IntNonlinMgr m_intNonlinMgr;
  AsymMgr      m_asymMgr;
  MPDMgr       m_mpdMgr;
  TholdCIMgr   m_tholdCIMgr;
  TholdMuonMgr m_tholdMuonMgr;
  
  /// hook the BeginEvent so that we can check our validity once per event.
  void handle ( const Incident& inc );

  //-- FRIEND CLASSES --//
  // following classes all share many properties w/ CalCalibSvc as they are
  // sort of 'employees' of CalCalibSvc.  easiest way to do it is to make them
  // friends
  friend class CalibItemMgr;
  friend class AsymMgr;
  friend class PedMgr;
  friend class IntNonlinMgr;
  friend class MPDMgr;
  friend class TholdCIMgr;
  friend class TholdMuonMgr;
};

#endif // CalCalibSvc_H
