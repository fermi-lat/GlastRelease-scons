#ifndef CalCalibSvc_H
#define CalCalibSvc_H
// $Header$
/** @file
    @author Z.Fewtrell
*/

// LOCAL 
#include "AsymMgr.h"
#include "IntNonlinMgr.h"
#include "MPDMgr.h"
#include "PedMgr.h"
#include "TholdCIMgr.h"
#include "CalCalibShared.h"

// GLAST 
#include "CalXtalResponse/ICalCalibSvc.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"

// EXTLIB
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IIncidentListener.h"

// STD

/** @class CalCalibSvc
    @author Z.Fewtrell
    * \brief Instatiates ICalCalibSvc interface, gets data from CalibDataSvc
    *
    * handles:
    * - data storage/destruction
    * - communication with Gleam lower level services
    * - checking of data validity period  
    * - extraction of cal-specific constants out of generic data objects
    * - creation/caching of local meta-data objects where needed (e.g. splines)
    *
    * jobOptions:
    * - CalibDataSvc - calibration data source (default="CalibDataSvc")
    * - idealCalibXMLPath - source for ideal calibration data (which by-passes the CalibSvc) (default="$(CALUTILXMLPATH)/idealCalib_flight.xml")
    * - DefaultFlavor - default flavor for all Cal Calibration types (default="ideal")
    * - FlavorIntNonlin - override default calibration flavor for IntNonlin only (default="" - disabled) 
    * - FlavorAsym - override default calib flavor for Asymmetry only (default="" - disabled)
    * - FlavorPed - override default calib flavor for pedestals only (default="" - disabled)
    * - FlavorMeVPerDac - override default calib flavor for mevPerDAC calibration (default="" - disabled)
    * - FlavorTholdCI - override default calib flavor for tholdCI calibration (default="" - disabled)
    *
    */

class CalCalibSvc : public Service, 
                    virtual public ICalCalibSvc,
                    virtual public IIncidentListener {
  
public:
  
  CalCalibSvc(const string& name, ISvcLocator* pSvcLocator); 
  
  StatusCode initialize();
  StatusCode finalize () {return StatusCode::SUCCESS;}

  /// queryInterface - for implementing a Service this is necessary
  StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);

  /// return the service type
  const InterfaceID&  type () const {return IID_ICalCalibSvc;}

  /// get calibration for given crystal
  const CalibData::CalMevPerDac *getMPD(CalUtil::XtalIdx xtalIdx) {return m_mpdMgr.getMPD(xtalIdx);}

  /// get adc points for given intNonlin (cidac2adc) calibration curve
  const std::vector<float> *getInlAdc(CalUtil::RngIdx rngIdx) {
    return m_inlMgr.getInlAdc(rngIdx);}

  /// get cidac points for given intNonlin calibration curve
  const std::vector<float> *getInlCIDAC(CalUtil::RngIdx rngIdx) {
    return m_inlMgr.getInlCIDAC(rngIdx);}

  /// get pedestal calibration data for given adc channel
  StatusCode getPed(CalUtil::RngIdx rngIdx, float &ped);
  StatusCode getPedSig(CalUtil::RngIdx rngIdx, float &sig);

  /// get asymmetry calib data for give crystal
  const CalibData::CalAsym *getAsym(CalUtil::XtalIdx xtalIdx) {return m_asymMgr.getAsym(xtalIdx);}
  
  /// get positions in mm from center of xtal for asymmetry
  /// calibratoin data points
  const CalibData::Xpos *getAsymXpos() {return m_asymMgr.getXpos();}

  /// get threshold calib values for given crystal face
  const CalibData::CalTholdCI *getTholdCI(CalUtil::FaceIdx faceIdx) {return m_tholdCIMgr.getTholdCI(faceIdx);}

  /// calculate associated cidac (charge-injection dac) signal level for given channel and adc value.
  StatusCode evalCIDAC (CalUtil::RngIdx rngIdx, float adc,   float &cidac) {
    return m_inlMgr.evalCIDAC(rngIdx, adc, cidac);
  }

  /// calculate associated adc value for given channel and cidac
  /// (charge-injection dac) setting
  StatusCode evalADC (CalUtil::RngIdx rngIdx, float cidac,   float &adc) {
    return m_inlMgr.evalADC(rngIdx, cidac, adc);
  }

  /// calculate light asymmetry value for given crystal and crystal position
  StatusCode evalAsym(CalUtil::XtalIdx xtalIdx, CalUtil::AsymType asymType, 
                      float pos,   float &asym) {
    return m_asymMgr.evalAsym(xtalIdx, asymType, pos, asym);
  }


  /// calculate crystal longitudinal position for energy deposit with
  /// given light asymmetry.
  StatusCode evalPos (CalUtil::XtalIdx xtalIdx, CalUtil::AsymType asymType, 
                      float asym,  float &pos) {
    return m_asymMgr.evalPos(xtalIdx, asymType, asym, pos);
  }

  /// calculate associated MeV deposit @ center of crystal that would
  /// result in given ADC readout for given ADC channel
  StatusCode evalFaceSignal(CalUtil::RngIdx rngIdx, float adc, float &ene);

  /// get MeVPerDAC value for single diode (instead of both faces)
  StatusCode getMPDDiode(CalUtil::DiodeIdx diodeIdx, float &mpdDiode);

  /// get light asymmetry associated with deposit at center of given crystal
  StatusCode getAsymCtr(CalUtil::XtalIdx xtalIdx, CalUtil::AsymType asymType, 
                        float &asymCtr) {
    return m_asymMgr.getAsymCtr(xtalIdx, asymType, asymCtr);
  }

  /// 'BeginEvent' handle needs lower priority (higher number) than CalibDataSvc
  static const unsigned short INCIDENT_PRIORITY = 50;

private:
  ////////////////////////////////////////////////
  ////// PARAMETER MANAGEMENT ////////////////////
  ////////////////////////////////////////////////

  // JobOptions PROPERTIES

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

  /// file with CU CAL tower temperature measurements
  StringProperty m_temperatureFile;
  /// file with pedestal temperature correction data
  StringProperty m_pedTempCorFile;

  ///CU temperature measurements table
  vector<int  > m_tempTime;
  vector<float> m_cuTwrTemp[4];
  float m_cur_temp_twr[4];
  int m_nEvent;

  ///CU pedestal temperature correction data
 
    CalUtil::CalVec<CalUtil::RngIdx, float> m_pedT0;
    CalUtil::CalVec<CalUtil::RngIdx, float> m_pedTempCoef;

  /// this class is shared amongt the CalibItemMgr classes
  CalCalibShared m_ccsShared;

  /// manage calibration data of given data type
  PedMgr       m_pedMgr;
  /// manage calibration data of given data type
  IntNonlinMgr m_inlMgr;
  /// manage calibration data of given data type
  AsymMgr      m_asymMgr;
  /// manage calibration data of given data type
  MPDMgr       m_mpdMgr;
  /// manage calibration data of given data type
  TholdCIMgr   m_tholdCIMgr;

  IDataProviderSvc* m_eventSvc;

  /// set to true when we have retrieved event time for current event.
  bool m_gotPedTime;


  /// hook the BeginEvent so that we can check our validity once per event.
  void handle ( const Incident& inc );

  /// get calib db serial number for current calibration data
  int getSerNoPed(){
    return m_pedMgr.getSerNo();
  }

  /// get calib db serial number for current calibration data
  int getSerNoINL(){
    return m_pedMgr.getSerNo();
  }

  /// get calib db serial number for current calibration data
  int getSerNoAsym(){
    return m_pedMgr.getSerNo();
  }

  /// get calib db serial number for current calibration data
  int getSerNoMPD(){
    return m_pedMgr.getSerNo();
  }

  /// get calib db serial number for current calibration data
  int getSerNoTholdCI(){
    return m_pedMgr.getSerNo();
  }
};

#endif // CalCalibSvc_H

