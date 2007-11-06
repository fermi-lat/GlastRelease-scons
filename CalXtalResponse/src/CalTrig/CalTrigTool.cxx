// $Header$

// Include files
/** @file
    @author Zach Fewtrell
*/

// LOCAL
#include "CalTrigTool.h"
#include "CalXtalResponse/ICalSignalTool.h"

// GLAST
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/Digi/GltDigi.h"
#include "Event/Digi/CalDigi.h"
#include "CalUtil/CalArray.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalGeom.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

// EXTLIB
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"

// STD
#include <algorithm>

using namespace CalUtil;
using namespace Event;
using namespace CalibData;
using namespace idents;


static ToolFactory<CalTrigTool> s_factory;
const IToolFactory& CalTrigToolFactory = s_factory;

CalTrigTool::CalTrigTool( const string& type,
                          const string& name,
                          const IInterface* parent)
  : AlgTool(type,name,parent),
    m_calCalibSvc(0),
    m_evtSvc(0),
    m_calSignalTool(0),
    m_detSvc(0)
{
  declareInterface<ICalTrigTool>(this);

  declareProperty("CalCalibSvc",      m_calCalibSvcName = "CalCalibSvc");
  declareProperty("PrecalcCalibTool", m_precalcCalibName = "PrecalcCalibTool");
  declareProperty("CalSignalToolName", m_calSignalToolName = "CalSignalTool");

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

  sc = toolSvc()->retrieveTool(m_calSignalToolName,
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
  m_twrList = CalUtil::findActiveTowers(*m_detSvc);

  //-- Retreive EventDataSvc
  sc = serviceLocator()->service( "EventDataSvc", m_evtSvc, true );
  if(sc.isFailure()){
    msglog << MSG::ERROR << "Could not find EventDataSvc" << endreq;
    return sc;
  }


  return StatusCode::SUCCESS;
}

/** 
    b/c we only have 1 range readout per face as input, the most consistent
    method to calculate both FLE & FHE triggers is by converting all adc
    readouts _and_ threshold levels to MeV at center of xtal & then we are always
    comparing the same units.  MeV readout -vs- MeV theshold regardless of 
    which adc range was recorded.
*/

StatusCode CalTrigTool::calcXtalTrig(XtalIdx xtalIdx,
                                     const CalDigi::CalXtalReadout &ro,
                                     CalArray<XtalDiode, bool> &trigBits,
                                     Event::GltDigi *glt
                                     ) {
  StatusCode sc;

  fill(trigBits.begin(), trigBits.end(), false);
  
  for (FaceNum face; face.isValid(); face++) {
    FaceIdx faceIdx(xtalIdx, face);
    XtalDiode xDiodeLrg(face, LRG_DIODE);
    XtalDiode xDiodeSm(face,  SM_DIODE);

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
    Ped const*const ped = m_calCalibSvc->getPed(rngIdx);
    if (!ped) return StatusCode::FAILURE;
    const float adcPed = adc - ped->getAvr();

    //-- eval faceSignal 
    sc = m_calCalibSvc->evalFaceSignal(rngIdx, adcPed, ene);
    if (sc.isFailure()) return sc;

    // set trigger bit
    if (ene >= fleMeV) {
      trigBits[xDiodeLrg] = true;
      // optionally populate GltDigi
      if (glt)
        glt->setCALLOtrigger(faceIdx.getCalXtalId());
    }

    // FHE //
    //-- RETRIEVE THESHOLD CALIB (per-face) --//
    float fheMeV;
    diodeIdx = DiodeIdx(faceIdx, SM_DIODE);
    sc = m_precalcCalibTool->getTrigMeV(diodeIdx, fheMeV);
    if (sc.isFailure()) return sc;

    // set trigger bit
    if (ene >= fheMeV) {
      trigBits[xDiodeSm] = true;
      // optionally populate GltDigi
      if (glt)
        glt->setCALHItrigger(faceIdx.getCalXtalId());
    }
  } // face loop
  
  return StatusCode::SUCCESS;
}


/**
   b/c we have all 4 adc range readouts, the stategy of this function is simple...
   for each threshold, simply check the appropriate adc readout for comparison.
*/
StatusCode CalTrigTool::calcXtalTrig(XtalIdx xtalIdx,
                                     const CalArray<XtalRng, float> &adcPed,
                                     CalArray<XtalDiode, bool> &trigBits,
                                     Event::GltDigi *glt
                                     ) {
  StatusCode sc;

  fill(trigBits.begin(), trigBits.end(), false);

  //-- RETRIEVE CALIB --//
  for (FaceNum face; face.isValid(); face++) {
    FaceIdx faceIdx(xtalIdx,face);
    XtalDiode xDiodeLrg(face, LRG_DIODE);
    XtalDiode xDiodeSm(face,  SM_DIODE);

    // FLE //
    DiodeIdx diodeIdx(faceIdx, LRG_DIODE);
    float fleADC;
    RngNum fleRng;
    sc = m_precalcCalibTool->getTrigRngADC(diodeIdx, fleRng, fleADC);
    if (sc.isFailure()) return sc;

    // set trigger bit
    if (adcPed[XtalRng(face,fleRng)] >= fleADC) {
      trigBits[xDiodeLrg] = true;
      // optionally populate GltDigi
      if (glt)
        glt->setCALLOtrigger(faceIdx.getCalXtalId());
    }

    // FHE //
    diodeIdx = DiodeIdx(faceIdx, SM_DIODE);
    float fheADC;
    RngNum fheRng;
    sc = m_precalcCalibTool->getTrigRngADC(diodeIdx, fheRng, fheADC);
    if (sc.isFailure()) return sc;

    // set trigger bit
    if (adcPed[XtalRng(face,fheRng)] >= fheADC) {
      trigBits[xDiodeSm] = true;
      // optionally populate GltDigi
      if (glt)
        glt->setCALHItrigger(faceIdx.getCalXtalId());
    }
  }

  return StatusCode::SUCCESS;
}

/**
   The purpose of this function is to call through to overloaded routines of the same 
   name which specialize in either 1, or 4 range readout for the GLAST Cal.
*/
StatusCode CalTrigTool::calcXtalTrig(const Event::CalDigi& calDigi,
                                     CalArray<XtalDiode, bool> &trigBits,
                                     Event::GltDigi *glt
                                     ) {
  XtalIdx xtalIdx(calDigi.getPackedId());

  //-- CASE 1: SINGLE READOUT MODE --//
  // if anything but 4 range mode, simply process 0th readout
  if (calDigi.getMode() != CalXtalId::ALLRANGE) {
    const Event::CalDigi::CalXtalReadout *ro = calDigi.getXtalReadout(0);
    if (!ro) {
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::ERROR << "Empty CalDigi (no '0' readout)" << endreq;
      return StatusCode::FAILURE;
    }

    return calcXtalTrig(xtalIdx, *ro, trigBits, glt);
  }

  //-- CASE 2: 4RANGE READOUT MODE --//
  else {
    //-- store ped subtracted adc vals --//
    CalArray<XtalRng, float> adcPed(0);

    //-- copy over CalDigi data from all available ranges --//
  
    for (FaceNum face; face.isValid(); face++)
      for (CalDigi::CalXtalReadoutCol::const_iterator ro = 
             calDigi.getReadoutCol().begin();
           ro != calDigi.getReadoutCol().end();
           ro++) {
        const RngNum rng((*ro).getRange(face));
      
        //-- retrieve pedestals --//
        RngIdx rngIdx(xtalIdx,
                      face, rng);

        Ped const*const ped = m_calCalibSvc->getPed(rngIdx);
        if (!ped) return StatusCode::FAILURE;
        
        adcPed[XtalRng(face, rng)]  = 
          (*ro).getAdc(face) - ped->getAvr();
      }

    return calcXtalTrig(xtalIdx, adcPed, trigBits, glt);
  }
}

StatusCode CalTrigTool::calcXtalTrig(const XtalIdx xtalIdx,
                                     const ICalSignalTool::XtalSignalMap &cidac,
                                     CalArray<XtalDiode, bool> &trigBits,
                                     Event::GltDigi *glt
                                     ) {
  StatusCode sc;

  fill(trigBits.begin(), trigBits.end(), false);

  // get threholds
  for (XtalDiode xDiode; xDiode.isValid(); xDiode++) {
    float thresh;
    sc = m_precalcCalibTool->getTrigCIDAC(DiodeIdx(xtalIdx, xDiode), thresh);
    if (sc.isFailure()) return sc;

    if (cidac[xDiode] >= thresh) {
      trigBits[xDiode] = true;
      
      // set glt digi info (optional)
      if (glt)
        if (xDiode.getDiode() == LRG_DIODE)
          glt->setCALLOtrigger(FaceIdx(xtalIdx, xDiode.getFace()).getCalXtalId());
        else
          glt->setCALHItrigger(FaceIdx(xtalIdx, xDiode.getFace()).getCalXtalId());
    }
  }
  
  return StatusCode::SUCCESS;
}
                                                                         
/** Either calc trigger from mc hits or rather from digi
*/
StatusCode CalTrigTool::calcGlobalTrig(CalArray<DiodeNum, bool> &trigBits,
                                       Event::GltDigi *glt) {
  // 1st check for presense of mc
  SmartDataPtr<Event::McIntegratingHitVector> mcHits(m_evtSvc, EventModel::MC::McIntegratingHitCol);
  if (mcHits != 0)
    return calcGlobalTrigSignalTool(trigBits,
                                    glt);
  
  // else calc with cal digis
  SmartDataPtr<Event::CalDigiCol> calDigis(m_evtSvc, EventModel::Digi::CalDigiCol);
  return calcGlobalTrig(calDigis, trigBits, glt);
}


/** Loop through digiCol & call calcXtalTrig() on each digi readout.
   Save to optional GltDigi if glt input param is non-null.
*/
StatusCode CalTrigTool::calcGlobalTrig(const Event::CalDigiCol &calDigiCol,
                                       CalArray<DiodeNum, bool> &trigBits,
                                       Event::GltDigi *glt) {
  StatusCode sc;

  //-- init --//
  fill(trigBits.begin(), trigBits.end(), false);

  CalArray<XtalDiode, bool> xtalTrigBits;
  
  // loop over all calorimeter digis in CalDigiCol
  for (CalDigiCol::const_iterator digiIter = calDigiCol.begin(); 
       digiIter != calDigiCol.end(); digiIter++) {
    sc = calcXtalTrig(**digiIter, xtalTrigBits, glt);
    if (sc.isFailure()) return sc;

    //-- set global trigger if any xtal trigger is high.
    for (XtalDiode xDiode; xDiode.isValid(); xDiode++) {
      const DiodeNum diode = xDiode.getDiode();
      
      trigBits[diode] |= xtalTrigBits[xDiode];
    }
  }

  return StatusCode::SUCCESS;
}


GltDigi* CalTrigTool::setupGltDigi() {
  StatusCode sc;

  // search for GltDigi in TDS
  static const string gltPath( EventModel::Digi::Event+"/GltDigi");
  DataObject* pnode = 0;
  auto_ptr<GltDigi> glt;

  sc = m_evtSvc->findObject(gltPath, pnode);
  // if the required entry doens't exit - create it
  if (sc.isFailure()) {
    glt.reset(new GltDigi());
    // always register glt data, even if there is no caldigi data.
    // sometimes we can trigger w/ no LACs.
    sc = m_evtSvc->registerObject(gltPath, glt.get());
    if (sc.isFailure()) {
      // if cannot create entry - error msg
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::WARNING << " Could not create GltDigi entry" << endreq;
      return NULL;
    }
  }
  else glt.reset(dynamic_cast<GltDigi*>(pnode));

  return glt.release();
}


/// loop through all crystals in installed towers, get signal level for each from CalSignalTool, calc trig resonse one by one
StatusCode CalTrigTool::calcGlobalTrigSignalTool(CalArray<DiodeNum, bool> &trigBits,
                                                 Event::GltDigi *glt) {
  fill(trigBits.begin(), trigBits.end(), false);

  /// Loop through (installed) towers and crystals; 
  for (unsigned twrSeq = 0; twrSeq < m_twrList.size(); twrSeq++) {
    // get bay id of nth live tower
    const TwrNum twr(m_twrList[twrSeq]);
    for (LyrNum lyr; lyr.isValid(); lyr++)
      for (ColNum col; col.isValid(); col++) {
        
        // assemble current calXtalId
        const XtalIdx xtalIdx(twr.val(),
                              lyr.val(),
                              col.val());
        ICalSignalTool::XtalSignalMap xtalSignalMap;
        StatusCode sc = m_calSignalTool->getXtalSignalMap(xtalIdx, xtalSignalMap);
        if (sc.isFailure())
          return sc;
        
        CalArray<XtalDiode, bool> xtalTrigBits;
        sc = calcXtalTrig(xtalIdx, xtalSignalMap, xtalTrigBits, glt);
        if (sc.isFailure())
          return sc;

        //-- set global trigger if any xtal trigger bit is high.
        for (XtalDiode xDiode; xDiode.isValid(); xDiode++) {
          const DiodeNum diode = xDiode.getDiode();
          
          trigBits[diode] |= xtalTrigBits[xDiode];
        }
      } // xtal loop
  } // twr loop

  return StatusCode::SUCCESS;
}


