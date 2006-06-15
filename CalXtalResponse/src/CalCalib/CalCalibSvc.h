#ifndef CalCalibSvc_H
#define CalCalibSvc_H
// $Header$

// LOCAL 
#include "AsymMgr.h"
#include "IntNonlinMgr.h"
#include "MPDMgr.h"
#include "PedMgr.h"
#include "TholdCIMgr.h"
#include "IdealCalCalib.h"
#include "CalCalibShared.h"

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
    virtual public IIncidentListener  {
  
    public:
  
  CalCalibSvc(const string& name, ISvcLocator* pSvcLocator); 
  
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize () {return StatusCode::SUCCESS;}



  /// queryInterface - for implementing a Service this is necessary
  StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);

  /// return the service type
  const InterfaceID&  CalCalibSvc::type () const {return IID_ICalCalibSvc;}

  const CalMevPerDac *getMPD(XtalIdx xtalIdx) {return m_mpdMgr.getMPD(xtalIdx);}

  const vector<float> *getInlAdc(CalUtil::RngIdx rngIdx) {
    return m_inlMgr.getInlAdc(rngIdx);}

  const vector<float> *getInlCIDAC(CalUtil::RngIdx rngIdx) {
    return m_inlMgr.getInlCIDAC(rngIdx);}

  const Ped *getPed(RngIdx rngIdx) {return m_pedMgr.getPed(rngIdx);}

  const CalAsym *getAsym(XtalIdx xtalIdx) {return m_asymMgr.getAsym(xtalIdx);}
  
  const Xpos *getAsymXpos() {return m_asymMgr.getXpos();}

  const CalTholdCI *getTholdCI(FaceIdx faceIdx) {return m_tholdCIMgr.getTholdCI(faceIdx);}

  StatusCode evalCIDAC (RngIdx rngIdx, float adc,   float &cidac) {
    return m_inlMgr.evalCIDAC(rngIdx, adc, cidac);
  }

  StatusCode evalADC (RngIdx rngIdx, float cidac,   float &adc) {
    return m_inlMgr.evalADC(rngIdx, cidac, adc);
  }

  StatusCode evalAsym(XtalIdx xtalIdx, AsymType asymType, 
                      float pos,   float &asym) {
    return m_asymMgr.evalAsym(xtalIdx, asymType, pos, asym);
  }

  StatusCode evalPos (XtalIdx xtalIdx, AsymType asymType, 
                      float asym,  float &pos) {
    return m_asymMgr.evalPos(xtalIdx, asymType, asym, pos);
  }

  StatusCode evalFaceSignal(RngIdx rngIdx, float adc, float &ene);

  StatusCode getAsymCtr(CalUtil::XtalIdx xtalIdx, CalUtil::AsymType asymType, 
                        float &asymCtr) {
    return m_asymMgr.getAsymCtr(xtalIdx, asymType, asymCtr);
  }


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


  /// this class is shared amongt the CalibItemMgr classes
  CalCalibShared m_ccsShared;

  PedMgr       m_pedMgr;
  IntNonlinMgr m_inlMgr;
  AsymMgr      m_asymMgr;
  MPDMgr       m_mpdMgr;
  TholdCIMgr   m_tholdCIMgr;

  /// hook the BeginEvent so that we can check our validity once per event.
  void handle ( const Incident& inc );

  int getSerNoPed(){
    return m_pedMgr.getSerNo();
  }

  int getSerNoINL(){
    return m_pedMgr.getSerNo();
  }

  int getSerNoAsym(){
    return m_pedMgr.getSerNo();
  }

  int getSerNoMPD(){
    return m_pedMgr.getSerNo();
  }

  int getSerNoTholdCI(){
    return m_pedMgr.getSerNo();
  }
};

#endif // CalCalibSvc_H
