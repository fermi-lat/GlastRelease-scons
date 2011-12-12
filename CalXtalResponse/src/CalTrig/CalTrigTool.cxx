// $Header$

// Include files

/** @file
    @author Z.Fewtrell
*/


// LOCAL
#include "CalTrigTool.h"
#include "CalXtalResponse/ICalSignalTool.h"
#include "../CalCalib/IPrecalcCalibTool.h"

// GLAST
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/Digi/CalDigi.h"
#include "CalUtil/CalVec.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalGeom.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

// EXTLIB
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentSvc.h"

// STD
#include <algorithm>

using namespace CalUtil;
using namespace Event;
using namespace CalibData;
using namespace idents;


//static ToolFactory<CalTrigTool> s_factory;
//const IToolFactory& CalTrigToolFactory = s_factory;
DECLARE_TOOL_FACTORY(CalTrigTool);

CalTrigTool::CalTrigTool( const std::string& type,
                          const std::string& name,
                          const IInterface* parent)
  : AlgTool(type,name,parent),
    m_calCalibSvc(0),
    m_precalcCalibTool(0),
    m_evtSvc(0),
    m_calSignalTool(0),
    m_detSvc(0),
    m_isValid(false),
    m_selectionRule("default")
{
  declareInterface<ICalTrigTool>(this);

  declareProperty("CalCalibSvc",       m_calCalibSvcName = "CalCalibSvc");
  declareProperty("PrecalcCalibTool",  m_precalcCalibName = "PrecalcCalibTool");
  declareProperty("CalSignalToolName", m_calSignalToolName = "CalSignalTool");
  declareProperty("selectionRule",     m_selectionRule= "default");
}

StatusCode CalTrigTool::initialize() {
  MsgStream msglog(msgSvc(), name());   
  msglog << MSG::INFO << "initialize" << endreq;

  StatusCode sc;

  //-- jobOptions --//
  if ((sc = setProperties()).isFailure()) {
    msglog << MSG::ERROR << "Failed to set properties" << endreq;
    return sc;
  }
    
  // obtain CalCalibSvc
  sc = service(m_calCalibSvcName.value(), m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get " << m_calCalibSvcName << endreq;
    return sc;
  }

  // this tool may also be shared by CalTrigTool, global ownership
  sc = toolSvc()->retrieveTool("PrecalcCalibTool", 
                               m_precalcCalibName, 
                               m_precalcCalibTool,
                               0); // shared
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  Unable to create " << m_precalcCalibName << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("CalSignalTool",
                               m_calSignalToolName,
                               m_calSignalTool,
                               0); // intended to be shared
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << m_calSignalToolName << endreq;
    return sc;
  }

  // now try to find the GlastDetSvc service
  sc = service("GlastDetSvc", m_detSvc);
  if (sc.isFailure() ) {
    MsgStream msglog(msgSvc(), name());
    msglog << MSG::ERROR << "  Unable to get GlastDetSvc " << endreq;
    return sc;
  }

  //-- find out which tems are installed.
  m_twrList = findActiveTowers(*m_detSvc);

  //-- Retreive EventDataSvc
  sc = serviceLocator()->service( "EventDataSvc", m_evtSvc, true );
  if(sc.isFailure()){
    msglog << MSG::ERROR << "Could not find EventDataSvc" << endreq;
    return sc;
  }

  // Get ready to listen for BeginEvent
  IIncidentSvc* incSvc;
  sc = service("IncidentSvc", incSvc, true);
  if (sc.isSuccess() ) {
    incSvc->addListener(this, "BeginEvent"); // priority not important
  } else {
    msglog << MSG::ERROR << "can't find IncidentSvc" << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}

StatusCode CalTrigTool::calcGlobalTrig() {
  // return success if data is valid for current event
  if (m_isValid)
    return StatusCode::SUCCESS;

  // 1st check for presense of mc
  SmartDataPtr<Event::McIntegratingHitVector> mcHits(m_evtSvc, EventModel::MC::McIntegratingHitCol);
  if (mcHits != 0) {
    if (calcGlobalTrigSignalTool().isFailure())
      return StatusCode::FAILURE;
  }
  
  // else calc with cal digis
  else {
    SmartDataPtr<Event::CalDigiCol> calDigiCol(m_evtSvc, EventModel::Digi::CalDigiCol);
    if (calDigiCol == 0) {
      MsgStream msglog(msgSvc(), name());   
      msglog << MSG::DEBUG << "Unable to retrieve cal digis (or MC) for cal trigger processing" << endreq;

      return StatusCode::SUCCESS;
    }

    if (calcGlobalTrigDigi(calDigiCol).isFailure())
      return StatusCode::FAILURE;
  }
  m_isValid = true;

  return StatusCode::SUCCESS;
}


/** Loop through digiCol & call calcXtalTrig() on each digi readout. */
StatusCode CalTrigTool::calcGlobalTrigDigi(const Event::CalDigiCol &calDigiCol) {
  StatusCode sc;

  // loop over all calorimeter digis in CalDigiCol
  for (CalDigiCol::const_iterator digiIter = calDigiCol.begin(); 
       digiIter != calDigiCol.end(); digiIter++) {
    sc = calcXtalTrig(**digiIter);
    if (sc.isFailure()) return sc;
  }

  return StatusCode::SUCCESS;
}


/// loop through all crystals in installed towers, get signal level for each from CalSignalTool, calc trig resonse one by one
StatusCode CalTrigTool::calcGlobalTrigSignalTool() {
  /// Loop through (installed) towers and crystals; 
  for (unsigned twrSeq = 0; twrSeq < m_twrList.size(); twrSeq++) {
    // get bay id of nth live tower
    const TwrNum twr(m_twrList[twrSeq]);
    for (LyrNum lyr; lyr.isValid(); lyr++)
      for (ColNum col; col.isValid(); col++) {
        
        // assemble current calXtalId
        const XtalIdx xtalIdx(twr,lyr, col);

        StatusCode sc = calcXtalTrigSignalTool(xtalIdx);
        if (sc.isFailure())
          return sc;
      } // xtal loop
  } // twr loop

  return StatusCode::SUCCESS;
}

StatusCode CalTrigTool::calcXtalTrigSignalTool(const XtalIdx xtalIdx) {
  StatusCode sc;

  // get threholds
  for (XtalDiode xDiode; xDiode.isValid(); xDiode++) {
    const DiodeIdx diodeIdx(xtalIdx, xDiode);

    /// get threshold
    float thresh;
    sc = m_precalcCalibTool->getTrigCIDAC(diodeIdx, thresh);
    if (sc.isFailure()) return sc;

    /// get signal level
    float signal;
    sc = m_calSignalTool->getDiodeSignal(diodeIdx, signal);
    if (sc.isFailure()) return sc;

    /// get signal level
    float trigSignal;
    sc = m_calSignalTool->getTrigDiodeSignal(diodeIdx, trigSignal);
    if (sc.isFailure()) return sc;

    if (trigSignal >= thresh)
      setSingleBit(diodeIdx);
  }
  
  return StatusCode::SUCCESS;
}

/**
   The purpose of this function is to call through to overloaded routines of the same 
   name which specialize in either 1, or 4 range readout for the GLAST Cal.
*/
StatusCode CalTrigTool::calcXtalTrig(const Event::CalDigi& calDigi) {
  const XtalIdx xtalIdx(calDigi.getPackedId());

  //-- CASE 1: SINGLE READOUT MODE --//
  // if anything but 4 range mode, simply process 0th readout
  if (calDigi.getMode() != CalXtalId::ALLRANGE) {
    Event::CalDigi::CalXtalReadout const*const ro = calDigi.getXtalReadout(0);
    if (!ro) {
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::ERROR << "Empty CalDigi (no '0' readout)" << endreq;
      return StatusCode::FAILURE;
    }

    return calcXtalTrig(xtalIdx, *ro);
  }

  //-- CASE 2: 4RANGE READOUT MODE --//
  else {
    //-- store ped subtracted adc vals --//
    CalVec<XtalRng, float> adcPed;

    //-- copy over CalDigi data from all available ranges --//
  
    for (FaceNum face; face.isValid(); face++)
      for (CalDigi::CalXtalReadoutCol::const_iterator ro = 
             calDigi.getReadoutCol().begin();
           ro != calDigi.getReadoutCol().end();
           ro++) {
        const RngNum rng((*ro).getRange(face));
      
        //-- retrieve pedestals --//
        const RngIdx rngIdx(xtalIdx,
                            face, rng);

    float ped;
    StatusCode sc = m_calCalibSvc->getPed(rngIdx,ped);
    if (sc.isFailure()) return StatusCode::FAILURE;

        adcPed[XtalRng(face, rng)]  = 
          (*ro).getAdc(face) - ped;
      }

    return calcXtalTrig(xtalIdx, adcPed);
  }
}

/// set all appropriate internal trigger bits
void CalTrigTool::setSingleBit(const CalUtil::DiodeIdx diodeIdx) {
    
    //----- check whether rows/columns selection requested from jobOptions
    // GCRCNum, ColNum
    
    short  crc= diodeIdx.getLyr().getGCRC().val();  // get "CAL row controller" number - can be 0-3 for a CAL face
    short  col= diodeIdx.getCol().val();            // get column number
    string sr=  m_selectionRule;                    // get rows/columns selection rule specified in jobOptions
        
    if(sr== "erec")
    {
        //----- disable trigger generation in [even rows, odd columns] OR [odd rows, even columns]
        
        if((crc%2== 0 && col%2!= 0) || (crc%2!= 0 && col%2== 0)) return;
    }
    else
    {
        if(sr== "eroc")
        {
            //----- disable trigger generation in [even rows, even columns] OR [odd rows, odd columns]
            
            if((crc%2== 0 && col%2== 0) || (crc%2!= 0 && col%2!= 0)) return;
        }
    }
    
    //----- full trigger map
            
    m_calTriggerMap[diodeIdx]= true;
        
    //----- cal trigger vectors
        
    const DiodeNum diode(diodeIdx.getDiode());
    const TwrNum   twr(diodeIdx.getTwr());
    m_calTriggerVec[diode] |= (1 << twr.val());
}


StatusCode CalTrigTool::getCALTriggerVector(const idents::CalXtalId::DiodeType diode, 
                                            unsigned short &vec) {
  /// update internal tables
  StatusCode sc(calcGlobalTrig());
  if (sc.isFailure())
    return sc;

  vec = m_calTriggerVec[DiodeNum(diode)];
  
  return StatusCode::SUCCESS;
}

StatusCode CalTrigTool::getTriggerBit(const CalUtil::DiodeIdx diodeIdx, bool &trigBit) {
  /// update internal tables
  StatusCode sc(calcGlobalTrig());
  if (sc.isFailure())
    return sc;

  trigBit = m_calTriggerMap[diodeIdx];

  return StatusCode::SUCCESS;
}


/** 
    b/c we only have 1 range readout per face as input, the most consistent
    method to calculate both FLE & FHE triggers is by converting all adc
    readouts _and_ threshold levels to MeV at center of xtal & then we are always
    comparing the same units.  MeV readout -vs- MeV theshold regardless of 
    which adc range was recorded.
*/
StatusCode CalTrigTool::calcXtalTrig(const XtalIdx xtalIdx,
                                     const CalDigi::CalXtalReadout &ro) {
  StatusCode sc;

  for (FaceNum face; face.isValid(); face++) {
    const FaceIdx faceIdx(xtalIdx, face);

    // FLE //
    //-- RETRIEVE THESHOLD CALIB (per-face) --//
    float fleMeV;
    DiodeIdx diodeIdx(faceIdx, LRG_DIODE);
    sc = m_precalcCalibTool->getTrigMeV(diodeIdx, fleMeV);
    if (sc.isFailure()) return sc;

    //-- CONVERT ADC READOUT TO MeV --//
    const RngNum rng(ro.getRange(face));
    const short adc = ro.getAdc(face);
    float ene;
    //-- retrieve pedestals 
    const RngIdx rngIdx(xtalIdx,
                        face, rng);

    float ped;
    StatusCode sc = m_calCalibSvc->getPed(rngIdx,ped);
    if (sc.isFailure()) return StatusCode::FAILURE;
    const float adcPed = adc - ped;

    //-- eval faceSignal 
    sc = m_calCalibSvc->evalFaceSignal(rngIdx, adcPed, ene);
    if (sc.isFailure()) return sc;

    // set trigger bit
    if (ene >= fleMeV)
      setSingleBit(diodeIdx);

    // FHE //
    //-- RETRIEVE THESHOLD CALIB (per-face) --//
    float fheMeV;
    diodeIdx = DiodeIdx(faceIdx, SM_DIODE);
    sc = m_precalcCalibTool->getTrigMeV(diodeIdx, fheMeV);
    if (sc.isFailure()) return sc;

    // set trigger bit
    if (ene >= fheMeV)
      setSingleBit(diodeIdx);
  } // face loop
  
  return StatusCode::SUCCESS;
}

/**
   b/c we have all 4 adc range readouts, the stategy of this function is simple...
   for each threshold, simply check the appropriate adc readout for comparison.
*/
StatusCode CalTrigTool::calcXtalTrig(const XtalIdx xtalIdx,
                                     const CalVec<XtalRng, float> &adcPed
                                     ) {
  StatusCode sc;

  //-- RETRIEVE CALIB --//
  for (FaceNum face; face.isValid(); face++) {
    const FaceIdx faceIdx(xtalIdx,face);

    // FLE //
    DiodeIdx diodeIdx(faceIdx, LRG_DIODE);
    float fleADC;
    RngNum fleRng;
    sc = m_precalcCalibTool->getTrigRngADC(diodeIdx, fleRng, fleADC);
    if (sc.isFailure()) return sc;

    // set trigger bit
    if (adcPed[XtalRng(face,fleRng)] >= fleADC)
      setSingleBit(diodeIdx);

    // FHE //
    diodeIdx = DiodeIdx(faceIdx, SM_DIODE);
    float fheADC;
    RngNum fheRng;
    sc = m_precalcCalibTool->getTrigRngADC(diodeIdx, fheRng, fheADC);
    if (sc.isFailure()) return sc;

    // set trigger bit
    if (adcPed[XtalRng(face,fheRng)] >= fheADC)
      setSingleBit(diodeIdx);
  }

  return StatusCode::SUCCESS;
}

/// Inform that a new incident has occured
void CalTrigTool::handle ( const Incident& inc ) { 
  if ((inc.type() == "BeginEvent"))
    newEvent();

  return; 
}

void CalTrigTool::newEvent() {
  fill(m_calTriggerMap.begin(),
       m_calTriggerMap.end(),
       false);

  fill(m_calTriggerVec.begin(),
       m_calTriggerVec.end(),
       0);

  m_isValid = false;
}

