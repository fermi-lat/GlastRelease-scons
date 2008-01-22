// $Header$
/** @file
    @author Z.Fewtrell
*/

// LOCAL INCLUDES
#include "test_CalDigi.h"
#include "CalXtalResponse/ICalCalibSvc.h"
#include "CalXtalResponse/ICalSignalTool.h"
#include "test_util.h"
#include "TestCfg.h"

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"

// EXTLIB INCLUDES

// STD INCLUDES
#include <set>

using namespace CalUtil;
using namespace std;
using namespace CalXtalResponse;

namespace {
  /// due to LACDAC pitch, atual threshold in mev may be incorrect.
  static const float LACDAC_MARGIN_MEV = .5;
}

/** 
    Algorithm:
    Allowing for some tolerance when determining zero suppression
    threshold response means that the exact hit-crystal list cannot
    not always be perfectly determined, so our tests become a little
    more vague.
    - check # of digis (within tolerance for case of zero suppression)
    - check that each xtal digi belongs to hit xtal in testCfg
    - check for duplicate xtals 
    (if total # of digis is correct & no duplicate xtals, & all digis
    belong to hit crystals, then we should be in good shape.

    - verify adc values for each hit xtal.
    
*/
StatusCode test_CalDigi::verify(ICalSignalTool &calSignalTool,
                                ICalCalibSvc &calCalibSvc,
                                const CalXtalResponse::TestCfg &testCfg,
                                const CalXtalResponse::TwrSet &twrSet,
                                const Event::CalDigiCol &calDigiCol) {
  

  //-- CHECK # OF DIGIS --//
  const unsigned short nDigis = calDigiCol.size();

  // CASE: zero supporess
  if (testCfg.zeroSuppress) {
    if (nDigis > testCfg.testXtals.size()) {
      MsgStream msglog(m_msgSvc, "test_CalDigi");   
      msglog << MSG::ERROR << "too many digis, zeroSuppress=true: " << endreq;
      return StatusCode::FAILURE;
    }

    const short xtalsOverLAC = testCfg.countXtalsOverLAC(calCalibSvc, twrSet, LACDAC_MARGIN_MEV);
    if (xtalsOverLAC < 0) //< failure case
      return StatusCode::FAILURE; 

    if (nDigis < xtalsOverLAC) {
      MsgStream msglog(m_msgSvc, "test_CalDigi");   
      msglog << MSG::ERROR << "too few digis, zeroSuppress=true: " << endreq;
      return StatusCode::FAILURE;
    }
  }

  // CASE: no zero suppress
  else {
    if (nDigis != twrSet.size()*tXtalIdx::N_VALS) {
      MsgStream msglog(m_msgSvc, "test_CalDigi");   
      msglog << MSG::ERROR << "wrong # of digis, zeroSuppress=false: " << endreq;
      return StatusCode::FAILURE;
    }
  }

  //-- check individual digis --//
  if (checkDuplicateXtals(calDigiCol).isFailure())
    return StatusCode::FAILURE;

  if (checkMissingXtals(calDigiCol, testCfg).isFailure())
    return StatusCode::FAILURE;

  // loop over all calorimeter digis in CalDigiCol
  for (Event::CalDigiCol::const_iterator digiIter = calDigiCol.begin(); 
       digiIter != calDigiCol.end(); digiIter++)
    if (verifyXtal(calSignalTool,
                   calCalibSvc,
                   testCfg,
                   **digiIter).isFailure())
      return StatusCode::FAILURE;


  

  return StatusCode::SUCCESS;
}

/// make sure each xtal digi belongs to a hit crystal in testCfg
StatusCode test_CalDigi::checkMissingXtals(const Event::CalDigiCol &calDigiCol,
                                           const CalXtalResponse::TestCfg &testCfg) {
  if (testCfg.zeroSuppress)
    for (Event::CalDigiCol::const_iterator digiIter = calDigiCol.begin(); 
         digiIter != calDigiCol.end(); digiIter++) {
      const XtalIdx xtalIdx((*digiIter)->getPackedId());
    
      if (find(testCfg.testXtals.begin(), testCfg.testXtals.end(), xtalIdx) == testCfg.testXtals.end()) {
        MsgStream msglog(m_msgSvc, "test_CalDigi");   
        msglog << MSG::ERROR << "wrong digi XtalId: " << xtalIdx.toStr() << endreq;
        return StatusCode::FAILURE;
      }
    }
  
  return StatusCode::SUCCESS;
}

StatusCode test_CalDigi::checkDuplicateXtals(const Event::CalDigiCol &calDigiCol) {
  CalVec<XtalIdx, bool> xtalSet;

  for (Event::CalDigiCol::const_iterator digiIter = calDigiCol.begin(); 
       digiIter != calDigiCol.end(); digiIter++) {
    const XtalIdx xtalIdx((*digiIter)->getPackedId());


    /// check that we have not already seen this xtal.
    if (xtalSet[xtalIdx]) {
      MsgStream msglog(m_msgSvc, "test_CalDigi");   
      msglog << MSG::ERROR << "duplicate xtal id in digi: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }
    xtalSet[xtalIdx] = true;
  }

  return StatusCode::SUCCESS;
  
}

/** Algorithm:
    - check correct number of readouts
    - check best range selection
    - check range ordering.
    - check ADC values
*/
StatusCode test_CalDigi::verifyXtal(ICalSignalTool &calSignalTool,
                                    ICalCalibSvc &calCalibSvc,
                                    const CalXtalResponse::TestCfg &testCfg,
                                    const Event::CalDigi &calDigi) {
  if (checkNReadouts(testCfg, calDigi).isFailure())
    return StatusCode::FAILURE;

  if (checkBestRange(calSignalTool, calCalibSvc, calDigi).isFailure())
    return StatusCode::FAILURE;

  if (checkRangeOrder(calDigi).isFailure())
    return StatusCode::FAILURE;

  if (checkADC(calSignalTool, calCalibSvc, calDigi).isFailure())
    return StatusCode::FAILURE;


  return StatusCode::SUCCESS;
}

StatusCode test_CalDigi::checkNReadouts(const CalXtalResponse::TestCfg &testCfg,
                                        const Event::CalDigi &calDigi) {
  unsigned short nExpectedReadouts;
  switch (testCfg.trigMode) {
  case idents::CalXtalId::ALLRANGE:
    nExpectedReadouts = 4;
    break;
  case idents::CalXtalId::BESTRANGE:
    nExpectedReadouts = 1;
    break;
  default:
    MsgStream msglog(m_msgSvc, "test_CalDigi");   
    msglog << MSG::ERROR << "Invalid Cal Trigger Mode: " << testCfg.trigMode << endreq;
    return StatusCode::FAILURE;
  } // swithc

  if (calDigi.getReadoutCol().size() != nExpectedReadouts) {
    MsgStream msglog(m_msgSvc, "test_CalDigi");   
    msglog << MSG::ERROR << "Wrong n readouts: " << calDigi.getReadoutCol().size() << endreq;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

StatusCode test_CalDigi::checkBestRange(ICalSignalTool &calSignalTool,
                                        ICalCalibSvc &calCalibSvc,
                                        const Event::CalDigi &calDigi
                                        ) {
  const XtalIdx xtalIdx(calDigi.getPackedId());

  // check each face independently
  for (FaceNum face; face.isValid(); face++) {
    const FaceIdx faceIdx(xtalIdx, face);

    /// find 1st unsaturated range for best range
    RngNum bestRng = LEX8; ///< will eventually hold 'best' range for this xtal face
    /// postcondition: if all ranges below HEX1 are saturated, bestRng will equal HEX1 when done.
    for (; bestRng < HEX1; bestRng++) {
      float dac = 0;
      const DiodeNum diode(bestRng.getDiode());
      const DiodeIdx diodeIdx(faceIdx, diode);
      if (calSignalTool.getDiodeSignal(diodeIdx, dac).isFailure())
        return StatusCode::FAILURE;

      float adc = 0;
      const RngIdx rngIdx(faceIdx, bestRng);
      if (calCalibSvc.evalADC(rngIdx, dac, adc).isFailure())
        return StatusCode::FAILURE;

      CalibData::CalTholdCI const*const tholdCI = calCalibSvc.getTholdCI(faceIdx);
      if (!tholdCI)
        return StatusCode::FAILURE;

      const float uldThresh = tholdCI->getULD(bestRng.val())->getVal();

      ///  if we are unsaturated, then we have found best range for
      ///  this xtal face
      if (adc < uldThresh)
        break;
    }

    /// check that digi matches best range
    if (calDigi.getRange(0, face) != bestRng.val()) {
      MsgStream msglog(m_msgSvc, "test_CalDigi");   
      msglog << MSG::ERROR << "Wrong best range: " << calDigi.getReadoutCol().size() << endreq;
      return StatusCode::FAILURE;
    }
    
  }

  return StatusCode::SUCCESS;
}

StatusCode test_CalDigi::checkRangeOrder(const Event::CalDigi &calDigi) {
  // only important in 4 range mode, (assume nReadouts is correct from
  // earlier tests)
  if (calDigi.getReadoutCol().size()  <= 1)
    return StatusCode::SUCCESS;

  /// 1st grab range info from best range
  CalArray<FaceNum, RngNum> lastRng;
  for (FaceNum face; face.isValid(); face++)
    lastRng[face] = calDigi.getRange(0,face);

  for (unsigned short nReadout = 1;
       nReadout < calDigi.getReadoutCol().size();
       nReadout++)
    for (FaceNum face; face.isValid(); face++) {
      const RngNum currentRng = calDigi.getRange(nReadout, face);
      if (currentRng.val() != (lastRng[face].val()+1) % RngNum::N_VALS) {
        MsgStream msglog(m_msgSvc, "test_CalDigi");   
        msglog << MSG::ERROR << "Wrong range order: " << endreq;
        return StatusCode::FAILURE;
      }
      lastRng[face] = currentRng;
    }

  return StatusCode::SUCCESS;
}

/** Algorithm:
    For each readout:
    - ped subtract the adc values from CalDigi
    - retrieve cidac values from calSignalTool
    - convert cidac signal to adc with calCalibSvc
    - clip testADC values to saturation thresholds
    - compare digi ADC values to test ADC values
*/
StatusCode test_CalDigi::checkADC(ICalSignalTool &calSignalTool,
                                  ICalCalibSvc &calCalibSvc,
                                  const Event::CalDigi &calDigi) {
  const XtalIdx xtalIdx(calDigi.getPackedId());

  for (Event::CalDigi::CalXtalReadoutCol::const_iterator ro =  calDigi.getReadoutCol().begin();
       ro != calDigi.getReadoutCol().end();
       ro++)
    for (FaceNum face; face.isValid(); face++) {
      const RngNum rng = ro->getRange(face);
      const float adc = ro->getAdc(face);
      const FaceIdx faceIdx(xtalIdx, face);

      // retrieve pedestal calibration
      const RngIdx rngIdx(faceIdx, rng);
      CalibData::Ped const*const pedCalib = calCalibSvc.getPed(rngIdx);
      if (!pedCalib) return 0;
      const float ped = pedCalib->getAvr();

      // ped subtracted ADC
      const float adcPed = adc - ped;

      /// get cidac value from calsignaltool
      float dac = 0;
      const DiodeNum diode(rng.getDiode());
      const DiodeIdx diodeIdx(faceIdx, diode);
      if (calSignalTool.getDiodeSignal(diodeIdx, dac).isFailure())
        return StatusCode::FAILURE;

      /// get adc from calcalibSvc
      float adcTest = 0;
      if (calCalibSvc.evalADC(rngIdx, dac, adcTest).isFailure())
        return StatusCode::FAILURE;

      // clip hex1 adc to uld value
      if (rng == HEX1) {
        CalibData::CalTholdCI const*const tholdCI = calCalibSvc.getTholdCI(faceIdx);
        if (!tholdCI)
          return StatusCode::FAILURE;

        const float uldThresh = tholdCI->getULD(rng.val())->getVal();
        adcTest = min(adcTest,uldThresh);
      }

      /// clip adcTest to MAX_ADC - ped
      adcTest = min(adcTest, 4095-ped);

      
      if (abs_diff(adcPed, adcTest) > MAX_ADC_DIFF) {
        MsgStream msglog(m_msgSvc, "test_CalDigi");   
        msglog << MSG::ERROR << "bad adc val: " << rngIdx.toStr() << endreq;
        return StatusCode::FAILURE;
      }
    }

  return StatusCode::SUCCESS;
}
