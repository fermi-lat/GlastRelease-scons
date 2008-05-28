// $Header$
/** @file 
    @author Z.Fewtrell
*/

// LOCAL INCLUDES
#include "test_CalCalibSvc.h"
#include "CalXtalResponse/ICalCalibSvc.h"
#include "test_util.h"
#include "TestCfg.h"

// GLAST INCLUDES
#include "CalUtil/stl_util.h"

// EXTLIB INCLUDES
#include "CLHEP/Random/RandFlat.h"

// STD INCLUDES
#include <string>
#include <cmath>

using namespace CalUtil;
using namespace std;
using namespace CalXtalResponse;

namespace {
  /// relative tolerance for calibration values (should only see
  /// rounding issues )
  static const float MAX_CALIB_DIFF = .001;

};

test_CalCalibSvc::TestCalibSet::TestCalibSet(const string &pedTXTPath,
                                             const string &cidac2adcTXTPath,
                                             const string &asymTXTPath,
                                             const string &mpdTXTPath,
                                             const string &tholdCITXTPath) 
{
  m_calPed.readTXT(pedTXTPath);
  m_cidac2adc.readTXT(cidac2adcTXTPath);
  m_calAsym.readTXT(asymTXTPath);
  m_calMPD.readTXT(mpdTXTPath);
  m_calTholdCI.readTXT(tholdCITXTPath);

}

/** Algorithm:
    - if appropriate tower module is missing, test for correct answer
    in this case
    - else test each major calibration type (ped, intNonlin,
    asymmetry, mevPerDAC, tholdCI)
*/
StatusCode test_CalCalibSvc::testXtal(const XtalIdx xtalIdx,
                                      ICalCalibSvc &calCalibSvc,
                                      const TestCalibSet &calibSet,
                                      const TwrSet &twrSet) {
  StatusCode sc;

  if (twrSet.find(xtalIdx.getTwr()) == twrSet.end())
    return testMissingXtal(xtalIdx,
                           calCalibSvc);
      
  sc = testPed(xtalIdx, calCalibSvc, calibSet);
  if (sc.isFailure())
    return sc;

  sc = testINL(xtalIdx, calCalibSvc, calibSet);
  if (sc.isFailure())
    return sc;

  sc = testAsym(xtalIdx, calCalibSvc, calibSet);
  if (sc.isFailure())
    return sc;

  sc = testMPD(xtalIdx, calCalibSvc, calibSet);
  if (sc.isFailure())
    return sc;

  sc = testTholdCI(xtalIdx, calCalibSvc, calibSet);
  if (sc.isFailure())
    return sc;

  if (testFaceSignal(xtalIdx, calCalibSvc).isFailure())
    return StatusCode::FAILURE;

  return StatusCode::SUCCESS;
}

/// check calCalibSvc values against TestCalibSet
StatusCode test_CalCalibSvc::testPed(const XtalIdx xtalIdx,
                                     ICalCalibSvc &calCalibSvc,
                                     const TestCalibSet &calibSet) {

  for (XtalRng xRng;
       xRng.isValid();
       xRng++) {
    const RngIdx rngIdx(xtalIdx, xRng);
    
    float val,sig;
    StatusCode sc = calCalibSvc.getPed(rngIdx,val);
    StatusCode scs = calCalibSvc.getPedSig(rngIdx,sig);
    if (sc.isFailure() || scs.isFailure()) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "missing pedestal: " << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }
    // check calibs against test values (should be close to w/in
    // rounding diff
    if (!smart_compare(val, calibSet.m_calPed.getPed(rngIdx), MAX_CALIB_DIFF)) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "Invalid pedestal: " << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    // check calibs against test values (should be close to w/in
    // rounding diff
    if (!smart_compare(sig, calibSet.m_calPed.getPedSig(rngIdx), MAX_CALIB_DIFF)) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "Invalid pedestal sigma: " << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }
 
  }

  return StatusCode::SUCCESS;
}

/** Algorithm:
   - check calCalibSvc values against TestCalibSet
   - check all points in spline (ADC & CIDAC)
   - for each point in spline, evaluate spline function both ways &
   check answer.
*/
StatusCode test_CalCalibSvc::testINL(const XtalIdx xtalIdx,
                                     ICalCalibSvc &calCalibSvc,
                                     const TestCalibSet &calibSet) {
  for (XtalRng xRng; xRng.isValid(); xRng++) {
    const RngIdx rngIdx(xtalIdx, xRng);

    vector<float> const*const inlCIDAC = calCalibSvc.getInlCIDAC(rngIdx);
    if (!inlCIDAC) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "missing INL CIDAC data: " << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    const vector<float> &cidacTest = calibSet.m_cidac2adc.getPtsDAC(rngIdx);
    if (!compare_collection(cidacTest.begin(),
                            cidacTest.end(),
                            inlCIDAC->begin(),
                            MAX_CALIB_DIFF)) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "invalid INL CIDAC data: " << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    vector<float> const*const inlADC = calCalibSvc.getInlAdc(rngIdx);
    if (!inlADC) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "missing INL ADC data: " << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    const vector<float> &adcTest = calibSet.m_cidac2adc.getPtsADC(rngIdx);
    if (!compare_collection(adcTest.begin(),
                            adcTest.end(),
                            inlADC->begin(),
                            MAX_CALIB_DIFF)) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "invalid INL ADC data: " << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    /// eval splines @ each point
    StatusCode sc;
    float tmp = 0;
    for (unsigned short i = 0;
         i < inlADC->size();
         i++) {

      sc = calCalibSvc.evalCIDAC(rngIdx,
                                 inlADC->at(i),
                                 tmp);
      if (sc.isFailure()) {
        MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
        msglog << MSG::ERROR << "INL spline failure: " << rngIdx.toStr() << endreq;
        return StatusCode::FAILURE;
      }

      if (!smart_compare(tmp, inlCIDAC->at(i),  MAX_SPLINE_DIFF)) {
        MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
        msglog << MSG::ERROR << "bad INL spline value: " << rngIdx.toStr() << endreq;
        return StatusCode::FAILURE;
      }

      sc = calCalibSvc.evalADC(rngIdx,
                               inlCIDAC->at(i),
                               tmp);
      if (sc.isFailure()) {
        MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
        msglog << MSG::ERROR << "INL spline failure: " << rngIdx.toStr() << endreq;
        return StatusCode::FAILURE;
      }

      if (!smart_compare(tmp, inlADC->at(i),  MAX_SPLINE_DIFF)) {
        MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
        msglog << MSG::ERROR << "bad INL spline value: " << rngIdx.toStr() << endreq;
        return StatusCode::FAILURE;
      }

      
    }

    
  }
  
  return StatusCode::SUCCESS;
}

/** check calCalibSvc values against TestCalibSet
    Algorithm:
    - check asymmetry value & xpos value for each point in each asym
    spline.
    - @ each point evalute spline function both ways & compare
    results against known answers.
    - test asymCtr method against evalAsym(pos=0)
 */
StatusCode test_CalCalibSvc::testAsym(const XtalIdx xtalIdx,
                                      ICalCalibSvc &calCalibSvc,
                                      const TestCalibSet &calibSet) {

  StatusCode sc;

  CalibData::CalAsym const*const asymCalib = calCalibSvc.getAsym(xtalIdx);
  if (!asymCalib) {
    MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
    msglog << MSG::ERROR << "missing Asym data: " << xtalIdx.toStr() << endreq;
    return StatusCode::FAILURE;
  }

  CalVec<AsymType, const vector<CalibData::ValSig> *> db_asym;
  db_asym[ASYM_LL] = asymCalib->getBig();
  db_asym[ASYM_SS] = asymCalib->getSmall();
  db_asym[ASYM_LS] = asymCalib->getNSmallPBig();
  db_asym[ASYM_SL] = asymCalib->getPSmallNBig();

  CalibData::Xpos const *const xpos = calCalibSvc.getAsymXpos();
  if (!xpos) {
    MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
    msglog << MSG::ERROR << "missing Asym xpos data: " << xtalIdx.toStr() << endreq;
    return StatusCode::FAILURE; 
  }
  const vector<float> *xvals = xpos->getVals();
  if (!xvals) {
    MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
    msglog << MSG::ERROR << "missing Asym xpos data: " << xtalIdx.toStr() << endreq;
    return StatusCode::FAILURE; 
  }

  for (AsymType asymType; asymType.isValid(); asymType++) {
    const vector<float> &asymTest = calibSet.m_calAsym.getPtsAsym(xtalIdx, asymType);
    const vector<float> &errTest  = calibSet.m_calAsym.getPtsErr(xtalIdx, asymType);

    if (db_asym[asymType]->size() == 0) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "missing Asym data: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    if (db_asym[asymType]->size() != asymTest.size()) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "invalid Asym data: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    for (unsigned short i = 0;
         i < asymTest.size();
         i++) {
      if (!smart_compare(db_asym[asymType]->at(i).getVal(), asymTest[i],  MAX_CALIB_DIFF)) {
        MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
        msglog << MSG::ERROR << "invalid Asym data: " << xtalIdx.toStr() << " " << asymType.toStr() << endreq;
        return StatusCode::FAILURE;
      }

      if (!smart_compare(db_asym[asymType]->at(i).getSig(), errTest[i],  MAX_CALIB_DIFF)) {
        MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
        msglog << MSG::ERROR << "invalid Asym data: " << xtalIdx.toStr() << " " << asymType.toStr() << endreq;
        return StatusCode::FAILURE;
      }

      //-- test spline methods
      float tmp = 0;
      sc = calCalibSvc.evalPos(xtalIdx, asymType, db_asym[asymType]->at(i).getVal(), tmp);
      if (sc.isFailure()) {
        MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
        msglog << MSG::ERROR << "Asym spline failure: " << xtalIdx.toStr() << endreq;
        return StatusCode::FAILURE;
      }
      
      if (!smart_compare(tmp, xvals->at(i), MAX_SPLINE_DIFF)) {
        MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
        msglog << MSG::ERROR << "bad Asym spline value: " << xtalIdx.toStr() << endreq;
        return StatusCode::FAILURE;
      }

      sc = calCalibSvc.evalAsym(xtalIdx, asymType, xvals->at(i), tmp);
      if (sc.isFailure()) {
        MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
        msglog << MSG::ERROR << "Asym spline failure: " << xtalIdx.toStr() << endreq;
        return StatusCode::FAILURE;
      }
      
      if (!smart_compare(tmp, db_asym[asymType]->at(i).getVal(), MAX_SPLINE_DIFF)) {
        MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
        msglog << MSG::ERROR << "bad Asym spline value: " << xtalIdx.toStr() << endreq;
        return StatusCode::FAILURE;
      }
      
    } // asym point loop

    //-- test AsymCtr method --//
    float tmp, tmp2 = 0;
    sc = calCalibSvc.getAsymCtr(xtalIdx,asymType, tmp);
    if (sc.isFailure()) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "getAsymCtr() failure: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    sc = calCalibSvc.evalAsym(xtalIdx, asymType, 0, tmp2);
    if (sc.isFailure()) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "Asym spline failure: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }
    
    if (!smart_compare(tmp, tmp2, MAX_SPLINE_DIFF)) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "getAsymCtr() failure: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }
    
  }

  return StatusCode::SUCCESS;
}

/// test calCalibSvc values against TestCalibSet
StatusCode test_CalCalibSvc::testMPD(const XtalIdx xtalIdx,
                                     ICalCalibSvc &calCalibSvc,
                                     const TestCalibSet &calibSet) {
  CalibData::CalMevPerDac const*const mpdCalib = calCalibSvc.getMPD(xtalIdx);
  if (!mpdCalib) {
    MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
    msglog << MSG::ERROR << "MISSING mevPerDAC: " << xtalIdx.toStr() << endreq;
    return StatusCode::FAILURE;
  }

  CalVec<DiodeNum, float> mpd, sig;
  mpd[CalUtil::LRG_DIODE] = mpdCalib->getBig()->getVal();
  mpd[CalUtil::SM_DIODE]  = mpdCalib->getSmall()->getVal();
  sig[CalUtil::LRG_DIODE] = mpdCalib->getBig()->getSig();
  sig[CalUtil::SM_DIODE]  = mpdCalib->getSmall()->getSig();

  for (DiodeNum diode; diode.isValid(); diode++) {
    if (!smart_compare(mpd[diode], calibSet.m_calMPD.getMPD(xtalIdx, diode),  MAX_CALIB_DIFF)) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "Invalid mevPerDAC: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    if (!smart_compare(sig[diode], calibSet.m_calMPD.getMPDErr(xtalIdx, diode),  MAX_CALIB_DIFF)) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "Invalid mevPerDAC sigma: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }
  }

  return StatusCode::SUCCESS;
}

/// test calCalibSvc values against TestCalibSet
StatusCode test_CalCalibSvc::testTholdCI(const XtalIdx xtalIdx,
                                         ICalCalibSvc &calCalibSvc,
                                         const TestCalibSet &calibSet) {

  for (FaceNum face; face.isValid(); face++) {
    const FaceIdx faceIdx(xtalIdx, face);

    CalibData::CalTholdCI const*const tholdCICalib = calCalibSvc.getTholdCI(faceIdx);
    if (!tholdCICalib) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "MISSING tholdCI: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    if (!smart_compare(tholdCICalib->getLAC()->getVal(), calibSet.m_calTholdCI.getLACThresh(faceIdx),  MAX_CALIB_DIFF)) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "Invalid LAC threshold: " << faceIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    if (!smart_compare(tholdCICalib->getFLE()->getVal(), calibSet.m_calTholdCI.getFLEThresh(faceIdx),  MAX_CALIB_DIFF)) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "Invalid FLE threshold: " << faceIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    if (!smart_compare(tholdCICalib->getFHE()->getVal(), calibSet.m_calTholdCI.getFHEThresh(faceIdx),  MAX_CALIB_DIFF)) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "Invalid FHE threshold: " << faceIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    for (RngNum rng; rng.isValid(); rng++) {
      const RngIdx rngIdx(faceIdx, rng);

      if (!smart_compare(tholdCICalib->getULD(rng.val())->getVal(), calibSet.m_calTholdCI.getULDThresh(rngIdx),  MAX_CALIB_DIFF)) {
        MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
        msglog << MSG::ERROR << "Invalid ULD threshold: " << faceIdx.toStr() << endreq;
        return StatusCode::FAILURE;
      }
    }

  }
  
  return StatusCode::SUCCESS;
}

/// ensure that correct answer is always returned in case of empty
/// tower bay for all methods in CalCalibSvc.
StatusCode test_CalCalibSvc::testMissingXtal(const XtalIdx xtalIdx,
                                             ICalCalibSvc &calCalibSvc) {
  StatusCode sc;

  //-- PER XTAL CALIBRATIONS --//
  if (calCalibSvc.getMPD(xtalIdx) != 0) {
    MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
    msglog << MSG::ERROR << "MPD calibrations for empty tower returned: " << xtalIdx.toStr() << endreq;
    return StatusCode::FAILURE;
  }

  if (calCalibSvc.getAsym(xtalIdx) != 0) {
    MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
    msglog << MSG::ERROR << "ASYM calibrations for empty tower returned: " << xtalIdx.toStr() << endreq;
    return StatusCode::FAILURE;
  }
   
  float tmp;
  //-- PER ASYM TYPE --//
  for (AsymType asymType;
       asymType.isValid();
       asymType++) {
    sc = calCalibSvc.evalAsym(xtalIdx, asymType, 0, tmp);
    if (sc.isSuccess()) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "ASYM calibrations for empty tower returned: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    sc = calCalibSvc.evalPos(xtalIdx, asymType, 0, tmp);
    if (sc.isSuccess()) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "ASYM calibrations for empty tower returned: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    sc = calCalibSvc.getAsymCtr(xtalIdx, asymType, tmp);
    if (sc.isSuccess()) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "ASYM calibrations for empty tower returned: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

  } // asym type

    //-- PER FACE CALIBRATIONS --//
  for (FaceNum face;
       face.isValid();
       face++) {
    const FaceIdx faceIdx(xtalIdx,face);

    if (calCalibSvc.getTholdCI(faceIdx) != 0) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "tholdCI calibrations for empty tower returned: " << faceIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }
     
  }

  //-- PER RANGE CALIBRATIONS --//
  for (XtalRng xRng;
       xRng.isValid();
       xRng++) {
    const RngIdx rngIdx(xtalIdx, xRng);
    float ped;
    if (!calCalibSvc.getPed(rngIdx,ped).isFailure()) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "Ped calibrations for empty tower returned: " << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    if (calCalibSvc.getInlAdc(rngIdx) != 0) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "INL calibrations for empty tower returned: " << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    if (calCalibSvc.getInlCIDAC(rngIdx) != 0) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "INL calibrations for empty tower returned: " << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    sc = calCalibSvc.evalCIDAC(rngIdx, 0, tmp);
    if (sc.isSuccess()) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "INL calibrations for empty tower returned: " << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    sc = calCalibSvc.evalADC(rngIdx, 0, tmp);
    if (sc.isSuccess()) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "INL calibrations for empty tower returned: " << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    sc = calCalibSvc.evalFaceSignal(rngIdx, 0, tmp);
    if (sc.isSuccess()) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "FaceSignal calibrations for empty tower returned: " << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }
  } // rngIdx loop
  return StatusCode::SUCCESS;
}


/// loop through each selected xtal in TestCfg & verify individutally.
StatusCode test_CalCalibSvc::verify(ICalCalibSvc &calCalibSvc,
                                    const test_CalCalibSvc::TestCalibSet &calibSet,
                                    const CalXtalResponse::TestCfg &testCfg,
                                    const TwrSet &twrSet) {
  StatusCode sc;
  MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   

  // xtal loop
  for (TestCfg::XtalList::const_iterator xtalIt(testCfg.testXtals.begin());
       xtalIt != testCfg.testXtals.end();
       xtalIt++) {

    sc = testXtal(*xtalIt, calCalibSvc, calibSet, twrSet);
    if (sc.isFailure())
      return sc;
  } // xtal loop
  return StatusCode::SUCCESS;
}
                                                   

/// test calCalibSvc->evalFaceSignal() against my own hand-calculated version.
StatusCode test_CalCalibSvc::testFaceSignal(const CalUtil::XtalIdx xtalIdx,
                                            ICalCalibSvc &calCalibSvc) {
  
  /// will need mevPerDAC for calculation
  CalibData::CalMevPerDac const*const mpdCalib = calCalibSvc.getMPD(xtalIdx);
  if (!mpdCalib) {
    MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
    msglog << MSG::ERROR << "MISSING mevPerDAC: " << xtalIdx.toStr() << endreq;
    return StatusCode::FAILURE;
  }

  CalVec<DiodeNum, float> mpd, sig;
  mpd[CalUtil::LRG_DIODE] = mpdCalib->getBig()->getVal();
  mpd[CalUtil::SM_DIODE]  = mpdCalib->getSmall()->getVal();

  /// test all adc channels in xtal
  for (XtalRng xRng;
       xRng.isValid();
       xRng++) {
    const RngIdx rngIdx(xtalIdx, xRng);
    
    /// get pedestal
    float ped;
    StatusCode sc = calCalibSvc.getPed(rngIdx,ped);
    if (sc.isFailure()) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "missing pedestal: " << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    /// select adc value within adc range.
    const float testADC = CLHEP::RandFlat::shoot(4095-ped);
    float testCIDAC(0);
    if (calCalibSvc.evalCIDAC(rngIdx, testADC, testCIDAC).isFailure())
      return StatusCode::FAILURE;
    

    /// get value to test
    float faceSignal(0);
    if (calCalibSvc.evalFaceSignal(rngIdx, testADC, faceSignal).isFailure())
      return StatusCode::FAILURE;

    // log(posDAC/negDAC)
    float asymCtr(0);
    const DiodeNum diode(xRng.getRng().getDiode());
    if (calCalibSvc.getAsymCtr(xtalIdx, 
                               AsymType(diode,diode),
                               asymCtr).isFailure())
      return StatusCode::FAILURE;

    // posDAC/negDAC
    const float asymRatio(exp(asymCtr));

    const FaceNum face(xRng.getFace());

    float meanCIDAC = 0;
    // asymRatio = POSDAC/NEGDAC
    // meanDAC = sqrt(POSDAC*NEGDAC)
    switch ((idents::CalXtalId::XtalFace)face) {
    case idents::CalXtalId::POS:
      // NEGDAC = POSDAC/asym
      // mean  = sqrt(POSDAC*POSDAC/asym)
      meanCIDAC = sqrt(testCIDAC*testCIDAC/asymRatio);
      break;
    case idents::CalXtalId::NEG:
      // POSDAC = NEGDAC*asym
      // mean = sqrt(NEGDAC*NEGDAC*asym)
      meanCIDAC = sqrt(testCIDAC*testCIDAC*asymRatio);
      break;
    default:
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "Invalid FaceNum" << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    const float testMeV = mpd[diode]*meanCIDAC;
    
    if (!smart_compare(testMeV, faceSignal, MAX_RECON_DIFF)) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "Invalid faceSignal: " << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

  }
       

  return StatusCode::SUCCESS;
}
