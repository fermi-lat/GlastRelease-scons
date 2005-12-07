// Include files

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


static ToolFactory<CalTrigTool> s_factory;
const IToolFactory& CalTrigToolFactory = s_factory;

CalTrigTool::CalTrigTool( const string& type,
                          const string& name,
                          const IInterface* parent)
  : AlgTool(type,name,parent),
    m_calCalibSvc(0)
{
  declareInterface<ICalTrigTool>(this);

  declareProperty("CalCalibSvc", m_calCalibSvcName = "CalCalibSvc");
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
    msglog << MSG::ERROR << "can't get CalCalibSvc." << endreq;
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

  m_dat.Clear();
  m_dat.twr = xtalIdx.getTwr();
  m_dat.lyr = xtalIdx.getLyr();
  m_dat.col = xtalIdx.getCol();

  trigBits.fill(false);
  
  for (FaceNum face; face.isValid(); face++) {
    FaceIdx faceIdx(xtalIdx, face);
    XtalDiode xDiodeLrg(face, LRG_DIODE);
    XtalDiode xDiodeSm(face,  SM_DIODE);

    //-- RETRIEVE THESHOLD CALIB (per-face) --//
    CalibData::ValSig fle, fhe, lac,
      uldTholdL8, uldTholdH8;
    sc = m_calCalibSvc->getTholdCI(faceIdx,fle,fhe,lac);
    if (sc.isFailure()) return sc;
    m_dat.trigThresh[xDiodeLrg] = fle.getVal();
    m_dat.trigThresh[xDiodeSm]  = fhe.getVal();

    sc = m_calCalibSvc->getULDCI(RngIdx(faceIdx,LEX8), uldTholdL8);
    if (sc.isFailure()) return sc;
    sc = m_calCalibSvc->getULDCI(RngIdx(faceIdx, HEX8), uldTholdH8);
    if (sc.isFailure()) return sc;

    m_dat.uldThold[XtalRng(face,LEX8)] = uldTholdL8.getVal();
    m_dat.uldThold[XtalRng(face,HEX8)] = uldTholdH8.getVal();


    //-- CONVERT ADC READOUT TO MeV --//
    RngNum rng(ro.getRange(face));
    short adc = ro.getAdc(face);
    float ene;

    //-- retrieve pedestals --//
    float ped;
    float sig, cos;  // placeholders
    RngIdx rngIdx(xtalIdx,
                  face, rng);

    sc = m_calCalibSvc->getPed(rngIdx,
                               ped,
                               sig, cos);
    if (sc.isFailure()) return sc;
    float adcPed = adc - ped;


    sc = m_calCalibSvc->evalFaceSignal(rngIdx, adcPed, ene);
    if (sc.isFailure()) return StatusCode::FAILURE;

    //-- FLE --//
    // thresholds are in X8 units, but may need conversion to X1 units if they
    // are past hardware range.
    RngNum fleRng = (m_dat.trigThresh[xDiodeLrg] > uldTholdL8.getVal()) ?
      LEX1 : LEX8;
    
    // convert to LEX1 range if needed
    if (fleRng == LEX1) {
      sc = lex8_to_lex1(faceIdx, 
                        m_dat.trigThresh[xDiodeLrg], 
                        m_dat.trigThresh[xDiodeLrg]);
      if (sc.isFailure()) return sc;
    }
    // convert FLE thresh to MeV
    float fleMeV;
    sc = m_calCalibSvc->evalFaceSignal(RngIdx(faceIdx,fleRng),
                                       m_dat.trigThresh[xDiodeLrg], fleMeV);
    if (sc.isFailure()) return sc;

    // set trigger bit
    if (ene >= fleMeV) {
      trigBits[xDiodeLrg] = true;
      // optionally populate GltDigi
      if (glt)
        glt->setCALLOtrigger(faceIdx.getCalXtalId());
    }


    //-- FHE --//
    // thresholds are in X8 units, but may need conversion to X1 units if they
    // are past hardware range.
    RngNum fheRng = (m_dat.trigThresh[xDiodeSm] > uldTholdH8.getVal()) ?
      HEX1 : HEX8;
    
    // convert to HEX1 range if needed
    if (fheRng == HEX1) {
      sc = hex8_to_hex1(faceIdx,                         
                        m_dat.trigThresh[xDiodeSm], 
                        m_dat.trigThresh[xDiodeSm]);

      if (sc.isFailure()) return sc;
    }
    // convert FHE thresh to MeV
    float fheMeV;
    sc = m_calCalibSvc->evalFaceSignal(RngIdx(faceIdx,fheRng),
                                       m_dat.trigThresh[xDiodeSm], fheMeV);

    // set trigger bit
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

  m_dat.Clear();
  m_dat.twr = xtalIdx.getTwr();
  m_dat.lyr = xtalIdx.getLyr();
  m_dat.col = xtalIdx.getCol();

  trigBits.fill(false);

  //-- RETRIEVE CALIB --//
  for (FaceNum face; face.isValid(); face++) {
    FaceIdx faceIdx(xtalIdx,face);
    XtalDiode xDiodeLrg(face, LRG_DIODE);
    XtalDiode xDiodeSm(face,  SM_DIODE);

    //-- THESHOLDS (per-face) --//
    CalibData::ValSig fle,fhe,lac,
      uldTholdL8, uldTholdH8;

    sc = m_calCalibSvc->getTholdCI(faceIdx,fle,fhe,lac);
    if (sc.isFailure()) return sc;
    m_dat.trigThresh[xDiodeLrg] = fle.getVal();
    m_dat.trigThresh[xDiodeSm]  = fhe.getVal();

    sc = m_calCalibSvc->getULDCI(RngIdx(faceIdx,LEX8),
                                 uldTholdL8);
    if (sc.isFailure()) return sc;
    sc = m_calCalibSvc->getULDCI(RngIdx(faceIdx,HEX8),
                                 uldTholdH8);
    if (sc.isFailure()) return sc;

    m_dat.uldThold[XtalRng(face,LEX8)] = uldTholdL8.getVal();
    m_dat.uldThold[XtalRng(face,HEX8)] = uldTholdH8.getVal();


    //-- FLE --//
    // thresholds are in X8 units, but may need conversion to X1 units if they
    // are past hardware range.
    RngNum fleRng = (m_dat.trigThresh[xDiodeLrg] > uldTholdL8.getVal()) ?
      LEX1 : LEX8;
    
    // convert to LEX1 range if needed
    if (fleRng == LEX1) {
      sc = lex8_to_lex1(faceIdx, 
                        m_dat.trigThresh[xDiodeLrg], 
                        m_dat.trigThresh[xDiodeLrg]);
      if (sc.isFailure()) return sc;
    }

    // set trigger bit
    if (adcPed[XtalRng(face,fleRng)] >= m_dat.trigThresh[xDiodeLrg] ) {
      trigBits[xDiodeLrg] = true;
      // optionally populate GltDigi
      if (glt)
        glt->setCALLOtrigger(faceIdx.getCalXtalId());
    }


    //-- FHE --//
    // thresholds are in X8 units, but may need conversion to X1 units if they
    // are past hardware range.
    RngNum fheRng = (m_dat.trigThresh[xDiodeSm] > uldTholdH8.getVal()) ?
      HEX1 : HEX8;
    
    // convert to HEX1 range if needed
    if (fheRng == HEX1) {
      sc = hex8_to_hex1(faceIdx,                         
                        m_dat.trigThresh[xDiodeSm], 
                        m_dat.trigThresh[xDiodeSm]);
      if (sc.isFailure()) return sc;
    }

    // set trigger bit
    if (adcPed[XtalRng(face,fheRng)] >= m_dat.trigThresh[xDiodeSm]) {
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
  StatusCode sc;

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
        float ped;
        float sig, cos;  // placeholders
        RngIdx rngIdx(xtalIdx,
                      face, rng);

        sc = m_calCalibSvc->getPed(rngIdx,
                                   ped,
                                   sig, cos);
        if (sc.isFailure()) return sc;
        
        adcPed[XtalRng(face, rng)]  = 
          (*ro).getAdc(face) - ped;
      }

    return calcXtalTrig(xtalIdx, adcPed, trigBits, glt);
  }
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
  GltDigi* glt=0;  // this will point to glt data one way or another
  DataObject* pnode = 0;
  sc = eventSvc->findObject(gltPath, pnode);
  glt = dynamic_cast<GltDigi*>(pnode);

  // if the required entry doens't exit - create it
  if (sc.isFailure()) {
    glt = new GltDigi();
    // always register glt data, even if there is no caldigi data.
    // sometimes we can trigger w/ no LACs.
    sc = eventSvc->registerObject(gltPath, glt);
    if (sc.isFailure()) {
      // if cannot create entry - error msg
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::WARNING << " Could not create GltDigi entry" << endreq;
      if (glt) delete glt;
      return NULL;
    }
  }

  return glt;
}

StatusCode CalTrigTool::lex8_to_lex1(FaceIdx faceIdx, float l8adc, float &l1adc) {

  float tmpDAC;
  StatusCode sc;

  FaceNum face = faceIdx.getFace();

  // use this point to generate LEX8/LEX1 ratio
  float x8tmp = m_dat.uldThold[XtalRng(face,LEX8)];

  // 1st convert to dac
  sc = m_calCalibSvc->evalDAC(RngIdx(faceIdx, LEX8),
                              x8tmp, tmpDAC);
  if (sc.isFailure()) return sc;
      
  // 2nd convert to LEX1 adc
  float x1tmp;
  sc = m_calCalibSvc->evalADC(RngIdx(faceIdx, LEX1),
                              tmpDAC, x1tmp);
  if (sc.isFailure()) return sc;
      
  float rat = x1tmp/x8tmp;

  l1adc = l8adc*rat;

  return StatusCode::SUCCESS;

}

StatusCode CalTrigTool::hex8_to_hex1(FaceIdx faceIdx, float h8adc, float &h1adc) {

  float tmpDAC;
  StatusCode sc;

  FaceNum face = faceIdx.getFace();

  // use this point to generate HEX8/HEX1 ratio
  float x8tmp = m_dat.uldThold[XtalRng(face,HEX8)];

  // 1st convert to dac
  sc = m_calCalibSvc->evalDAC(RngIdx(faceIdx, HEX8),
                              x8tmp, tmpDAC);
  if (sc.isFailure()) return sc;
      
  // 2nd convert to LEX1 adc
  float x1tmp;
  sc = m_calCalibSvc->evalADC(RngIdx(faceIdx, HEX1),
                              tmpDAC, x1tmp);
  if (sc.isFailure()) return sc;
      
  float rat = x1tmp/x8tmp;

  h1adc = h8adc*rat;

  return StatusCode::SUCCESS;

}

StatusCode CalTrigTool::finalize() {
  return StatusCode::SUCCESS;
}

