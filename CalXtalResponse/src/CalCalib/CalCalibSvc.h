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

using namespace CalUtil;

/** @class CalCalibSvc
    @author Zachary Fewtrell
    * \brief Instatiates ICalCalibSvc interface, gets data from CalibDataSvc
    *
    * handles:
    * - data storage/destruction
    * - communication with Gleam lower level services
    * - checking of data validity period  
    * - extraction of cal-specific constants out of generic data objects
    * - creation/caching of local meta-data objects where needed (e.g. splines)
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

  /// return the service type
  const IID&  type () const {return IID_ICalCalibSvc;}

  /// get MeVPerDac ratios for given xtal
  StatusCode getMPD(XtalIdx xtalIdx,
                    CalibData::ValSig &mpdLrg,
                    CalibData::ValSig &mpdSm) {
    return m_mpdMgr.getMPD(xtalIdx, mpdLrg, mpdSm);
  }

  /// get integral non-linearity vals for given xtal/face/rng
  StatusCode getIntNonlin(RngIdx rngIdx,
                          const vector< float > *&adcs,
                          const vector< float > *&dacs,
                          float &error) {
    return m_intNonlinMgr.getIntNonlin(rngIdx, adcs, dacs, error);
  }

  /// get pedestal vals for given xtal/face/rng
  StatusCode getPed(RngIdx rngIdx,
                    float &avr,
                    float &sig,
                    float &cos) {
    return m_pedMgr.getPed(rngIdx, avr, sig, cos);
  }

  /// get Asymmetry calibration information for one xtal
  StatusCode getAsym(XtalIdx xtalIdx,
                     const vector<CalibData::ValSig> *&asymLrg,
                     const vector<CalibData::ValSig> *&asymSm,
                     const vector<CalibData::ValSig> *&AsymNSPB,
                     const vector<CalibData::ValSig> *&asymPSNB,
                     const vector<float> *&xVals) {
    return m_asymMgr.getAsym(xtalIdx, asymLrg, asymSm, AsymNSPB, asymPSNB, xVals);
  }

  /// get threshold calibration constants as measured w/ charnge injection
  StatusCode getTholdCI(FaceIdx faceIdx,
                        CalibData::ValSig &FLE,
                        CalibData::ValSig &FHE,
                        CalibData::ValSig &LAC
                        ) {
    return m_tholdCIMgr.getTholds(faceIdx, FLE, FHE, LAC);
  }

  /// get Upper Level Discriminator threshold as measured w/ charnge
  /// injection for given xtal/face/rng

  StatusCode getULDCI(RngIdx rngIdx,
                      CalibData::ValSig &ULDThold) {
    return m_tholdCIMgr.getULD(rngIdx, ULDThold);
  }

  /// get pedestal calibration constants as measured during charge
  /// injection threshold testing.

  StatusCode getPedCI(RngIdx rngIdx,
                      CalibData::ValSig &ped) {
    return m_tholdCIMgr.getPed(rngIdx, ped);
  }

  /// get threshold calibration constants as measured w/ muon calibration
  StatusCode getTholdMuon(FaceIdx faceIdx,
                          CalibData::ValSig &FLE,
                          CalibData::ValSig &FHE
                          ) {
    return m_tholdMuonMgr.getTholds(faceIdx, FLE, FHE);
  }

  /// get pedestal calibration constants as measured during muon
  /// calibration threshold testing.

  StatusCode getPedMuon(RngIdx rngIdx,
                        CalibData::ValSig &ped) {
    return m_tholdMuonMgr.getPed(rngIdx, ped);
  }




  StatusCode getMPD(XtalIdx xtalIdx, 
                    CalArray<DiodeNum, float> &mpd);

  StatusCode getPed(XtalIdx xtalIdx,
                    CalArray<XtalRng, float> &peds,
                    CalArray<XtalRng, float> &sigs);

  StatusCode getTholdCI(XtalIdx xtalIdx,
                        CalArray<XtalDiode, float> &trigThesh,
                        CalArray<FaceNum, float> &lacThresh);

  StatusCode getULDCI(XtalIdx xtalIdx,
                      CalArray<XtalRng, float> &uldThold);



  StatusCode evalDAC(RngIdx rngIdx, float adc, float &dac) {
    return m_intNonlinMgr.evalDAC(rngIdx, adc, dac);
  }
  
  StatusCode evalADC(RngIdx rngIdx, float dac, float &adc) {
    return m_intNonlinMgr.evalADC(rngIdx, dac, adc);
  }
  
  StatusCode evalAsymLrg(XtalIdx xtalIdx, float pos, float &asymLrg) {
    return m_asymMgr.evalAsymLrg(xtalIdx, pos, asymLrg);
  }
  
  StatusCode evalPosLrg(XtalIdx xtalIdx, float asymLrg, float &pos) {
    return m_asymMgr.evalPosLrg(xtalIdx, asymLrg, pos);
  }
  
  StatusCode evalAsymSm(XtalIdx xtalIdx, float pos, float &asymSm) {
    return m_asymMgr.evalAsymSm(xtalIdx, pos, asymSm);
  }
  
  StatusCode evalPosSm(XtalIdx xtalIdx, float asymSm, float &pos) {
    return m_asymMgr.evalPosSm(xtalIdx, asymSm, pos);
  }
  
  StatusCode evalAsymNSPB(XtalIdx xtalIdx, float pos, float &asymNSPB) {
    return m_asymMgr.evalAsymNSPB(xtalIdx, pos, asymNSPB);
  }
  
  StatusCode evalPosNSPB(XtalIdx xtalIdx, float asymNSPB, float &pos) {
    return m_asymMgr.evalPosNSPB(xtalIdx, asymNSPB, pos);
  }
  
  StatusCode evalAsymPSNB(XtalIdx xtalIdx, float pos, float &asymPSNB) {
    return m_asymMgr.evalAsymPSNB(xtalIdx, pos, asymPSNB);
  }
  
  StatusCode evalPosPSNB(XtalIdx xtalIdx, float asymPSNB, float &pos) {
    return m_asymMgr.evalPosPSNB(xtalIdx, asymPSNB, pos);
  }

  StatusCode evalFaceSignal(RngIdx rngIdx, float adcPed, float &ene);


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

  // GAUDI RESOURCES
  /// pointer to CalibDataSvc
  IService         *m_calibDataSvc;     

  /// pointer to IDataProviderSvc interface of CalibDataSvc
  IDataProviderSvc *m_dataProviderSvc;  

  /// stores 'ideal' flavor generic calibration constants
  IdealCalCalib m_idealCalib;

  /// Load ideal calibration values from xml file
  StatusCode loadIdealCalib() {
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
