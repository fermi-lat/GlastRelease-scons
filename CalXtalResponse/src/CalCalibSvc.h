#ifndef CalCalibSvc_H
#define CalCalibSvc_H 1

// LOCAL
#include "AsymMgr.h"
#include "IntNonlinMgr.h"
#include "MPDMgr.h"
#include "PedMgr.h"
#include "TholdCIMgr.h"
#include "TholdMuonMgr.h"

// GLAST
#include "CalXtalResponse/ICalCalibSvc.h"

// EXTLIB
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

class CalCalibSvc : public Service, virtual public ICalCalibSvc, virtual public IIncidentListener {

 public:

  CalCalibSvc(const string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

  /// queryInterface - for implementing a Service this is necessary
  StatusCode queryInterface(const IID& riid, void** ppvUnknown);

  static const InterfaceID& interfaceID() {return ICalCalibSvc::interfaceID();}

  /// return the service type
  const IID& type() const;

  /// retrieve MeVPerDac ratios for given xtal
  StatusCode getMeVPerDac(const CalXtalId &xtalId,
                          CalibData::ValSig &lrg,
                          CalibData::ValSig &sm) {
    return m_mpdMgr.getMPD(xtalId, lrg, sm);
  }

  /// retrieve integral non-linearity values for given xtal/face/rng
  StatusCode getIntNonlin(const CalXtalId &xtalId,
                          const vector< float > *&adcs,
                          const vector< unsigned > *&dacs,
                          float &error) {
    return m_intNonlinMgr.getIntNonlin(xtalId, adcs, dacs, error);
  }

  /// retrieve pedestal values for given xtal/face/rng
  StatusCode getPed(const CalXtalId &xtalId,
                    float &avr,
                    float &sig,
                    float &cos) {
    return m_pedMgr.getPed(xtalId, avr, sig, cos);
  }

  /// retrieve Asymmetry calibration information for one xtal
  StatusCode getAsym(const CalXtalId &xtalId,
                     const vector<CalibData::ValSig> *&lrg,
                     const vector<CalibData::ValSig> *&sm,
                     const vector<CalibData::ValSig> *&NSmPLrg,
                     const vector<CalibData::ValSig> *&PSmNLrg,
                     const vector<float> *&xVals) {
    return m_asymMgr.getAsym(xtalId, lrg, sm, NSmPLrg, PSmNLrg, xVals);
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
  
  StatusCode evalAsymLrg(const CalXtalId &xtalId, double pos, double &lrg) {
    return m_asymMgr.evalAsymLrg(xtalId, pos, lrg);
  }
  
  StatusCode evalPosLrg(const CalXtalId &xtalId, double lrg, double &pos) {
    return m_asymMgr.evalPosLrg(xtalId, lrg, pos);
  }
  
  StatusCode evalAsymSm(const CalXtalId &xtalId, double pos, double &sm) {
    return m_asymMgr.evalAsymSm(xtalId, pos, sm);
  }
  
  StatusCode evalPosSm(const CalXtalId &xtalId, double sm, double &pos) {
    return m_asymMgr.evalPosSm(xtalId, sm, pos);
  }
  
  StatusCode evalAsymNSPB(const CalXtalId &xtalId, double pos, double &nspb) {
    return m_asymMgr.evalAsymNSPB(xtalId, pos, nspb);
  }
  
  StatusCode evalPosNSPB(const CalXtalId &xtalId, double nspb, double &pos) {
    return m_asymMgr.evalPosNSPB(xtalId, nspb, pos);
  }
  
  StatusCode evalAsymPSNB(const CalXtalId &xtalId, double pos, double &psnb) {
    return m_asymMgr.evalAsymPSNB(xtalId, pos, psnb);
  }
  
  StatusCode evalPosPSNB(const CalXtalId &xtalId, double psnb, double &pos) {
    return m_asymMgr.evalPosPSNB(xtalId, psnb, pos);
  }

 private:
  ////////////////////////////////////////////////
  ////// PARAMETER MANAGEMENT ////////////////////
  ////////////////////////////////////////////////

  // JobOptions PROPERTIES
  StringProperty m_calibDataSvcName;     ///< name of CalibDataSvc, main data source
  StringProperty m_defaultFlavor;        ///<  default flavor for all calib types, unless otherwise specified.

  // calib_item specific flavors override defaultFlavor
  StringProperty m_flavorIntNonlin;      ///< calib flavor override for int-nonlin constants
  StringProperty m_flavorAsym;           ///< calib flavor override for asymmetry constants
  StringProperty m_flavorPed;            ///< calib flavor override for pedestal constants
  StringProperty m_flavorMPD;            ///< calib flavor override for MeVPerDac constants
  StringProperty m_flavorTholdCI;        ///< calib flavor override for CI measured thresholds
  StringProperty m_flavorTholdMuon;      ///< calib flavor override for Muon measure thresholds

  // GAUDI RESOURCES
  IService         *m_calibDataSvc;     ///< pointer to CalibDataSvc
  IDataProviderSvc *m_dataProviderSvc;  ///< pointer to IDataProviderSvc interface of CalibDataSvc

  MPDMgr       m_mpdMgr;
  PedMgr       m_pedMgr;
  AsymMgr      m_asymMgr;
  IntNonlinMgr m_intNonlinMgr;
  TholdMuonMgr m_tholdMuonMgr;
  TholdCIMgr   m_tholdCIMgr;

  /// hook the BeginEvent so that we can check our validity once per event.
  void handle ( const Incident& inc );
};

#endif // CalCalibSvc_H

