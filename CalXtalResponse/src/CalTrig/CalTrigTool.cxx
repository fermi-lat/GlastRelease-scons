// $Header$

// Include files
/** @file
    @author Zach Fewtrell
 */

// LOCAL
#include "CalTrigTool.h"

// GLAST
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/EventModel.h"

// EXTLIB
#include "GaudiKernel/ToolFactory.h"

#include "GaudiKernel/MsgStream.h"

// STD

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
    m_calCalibSvc(0)
{
  declareInterface<ICalTrigTool>(this);

  declareProperty("CalCalibSvc",      m_calCalibSvcName = "CalCalibSvc");
  declareProperty("PrecalcCalibTool", m_precalcCalibName = "PrecalcCalibTool");
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
  sc = toolSvc()->retrieveTool("PrecalcCalibTool", m_precalcCalibName, m_precalcCalibTool);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  Unable to create " << m_precalcCalibName << endreq;
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

  trigBits.fill(false);
  
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
    RngNum rng(ro.getRange(face));
    short adc = ro.getAdc(face);
    float ene;
    //-- retrieve pedestals 
    RngIdx rngIdx(xtalIdx,
                  face, rng);
    const Ped *ped = m_calCalibSvc->getPed(rngIdx);
    if (!ped) return StatusCode::FAILURE;
    float adcPed = adc - ped->getAvr();
    //-- eval faceSignal 
    sc = m_calCalibSvc->evalFaceSignal(rngIdx, adcPed, ene);
    if (sc.isFailure()) return StatusCode::FAILURE;

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

  trigBits.fill(false);

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
   The purpose of this function is to callthrough to overloaded routines of the same 
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
    CalArray<XtalRng, float> adcPed;
    // "-1" indicates no data for given readout range
    adcPed.fill(0);

    //-- copy over CalDigi data from all available ranges --//
  
    for (FaceNum face; face.isValid(); face++)
      for (CalDigi::CalXtalReadoutCol::const_iterator ro = 
             calDigi.getReadoutCol().begin();
           ro != calDigi.getReadoutCol().end();
           ro++) {
        RngNum rng((*ro).getRange(face));
      
        //-- retrieve pedestals --//
        RngIdx rngIdx(xtalIdx,
                      face, rng);

        const Ped *ped = m_calCalibSvc->getPed(rngIdx);
        if (!ped) return StatusCode::FAILURE;
        
        adcPed[XtalRng(face, rng)]  = 
          (*ro).getAdc(face) - ped->getAvr();
      }

    return calcXtalTrig(xtalIdx, adcPed, trigBits, glt);
  }
}

StatusCode CalTrigTool::calcXtalTrig(XtalIdx xtalIdx,
                                     const CalArray<XtalDiode, float> &cidac,
                                     CalArray<XtalDiode, bool> &trigBits,
                                     Event::GltDigi *glt
                                     ) {
  StatusCode sc;

  trigBits.fill(false);

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
									 
/**
   Loop through digiCol & call calcXtalTrig() on each digi readout.
   Save to optional GltDigi if glt input param is non-null.

*/

StatusCode CalTrigTool::calcGlobalTrig(const CalDigiCol& calDigiCol,
                                       CalArray<DiodeNum, bool> &trigBits,
                                       Event::GltDigi *glt) {
  StatusCode sc;

  //-- init --//
  trigBits.fill(false);

  CalArray<XtalDiode, bool> xtalTrigBits;
  
  // loop over all calorimeter digis in CalDigiCol
  for (CalDigiCol::const_iterator digiIter = calDigiCol.begin(); 
       digiIter != calDigiCol.end(); digiIter++) {
    sc = calcXtalTrig(**digiIter, xtalTrigBits, glt);
    if (sc.isFailure()) return sc;

    //-- set global trigger if any xtal trigger is high.
    for (XtalDiode xDiode; xDiode.isValid(); xDiode++) {
      DiodeNum diode = xDiode.getDiode();
      
      trigBits[diode] |= xtalTrigBits[xDiode];
    }
  }

  return StatusCode::SUCCESS;
}


GltDigi* CalTrigTool::setupGltDigi(IDataProviderSvc *eventSvc) {
  StatusCode sc;

  // search for GltDigi in TDS
  static const string gltPath( EventModel::Digi::Event+"/GltDigi");
  DataObject* pnode = 0;
  sc = eventSvc->findObject(gltPath, pnode);

  auto_ptr<GltDigi> glt;
  // if the required entry doens't exit - create it
  if (sc.isFailure()) {
    glt.reset(new GltDigi());
    // always register glt data, even if there is no caldigi data.
    // sometimes we can trigger w/ no LACs.
    sc = eventSvc->registerObject(gltPath, glt.get());
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


