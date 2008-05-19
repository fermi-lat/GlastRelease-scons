// $Header$
/** @file 
    @author Z.Fewtrell
*/

// LOCAL INCLUDES
#include "test_PrecalcCalibTool.h"
#include "../CalCalib/IPrecalcCalibTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"
#include "test_util.h"
#include "TestCfg.h"

// GLAST INCLUDES


// EXTLIB INCLUDES
#include "GaudiKernel/MsgStream.h"

// STD INCLUDES
#include <string>
#include <cmath>

using namespace CalUtil;
using namespace std;
using namespace CalXtalResponse;


/// basically test each PrecalcCalibTool function against values we
/// calculate by hand from calCalibSvc
StatusCode test_PrecalcCalibTool::testXtal(const XtalIdx xtalIdx,
                                           ICalCalibSvc &calCalibSvc,
                                           IPrecalcCalibTool &precalcCalibTool) {

  //-- TEST getPedCIDAC() --//
  for (XtalRng xRng; xRng.isValid(); xRng++) {
    const RngIdx rngIdx(xtalIdx, xRng);
    float pedSigCIDAC;
    if (precalcCalibTool.getPedSigCIDAC(rngIdx, pedSigCIDAC).isFailure())
      return StatusCode::FAILURE;

    float pedSigADC;
    StatusCode sc = calCalibSvc.getPedSig(rngIdx,pedSigADC);
    if (sc.isFailure()) {
      MsgStream msglog(m_msgSvc, "test_PrecalcCalibTool");   
      msglog << MSG::ERROR << "missing pedestal: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    float testSigCIDAC;
    if (calCalibSvc.evalCIDAC(rngIdx, pedSigADC, testSigCIDAC).isFailure())
      return StatusCode::FAILURE;

    if (!smart_compare(pedSigCIDAC, testSigCIDAC, MAX_SPLINE_DIFF)) {
      MsgStream msglog(m_msgSvc, "test_PrecalcCalibTool");   
      msglog << MSG::ERROR << "Invalid pedestal sigma: " << rngIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }
  }

  //-- TEST getTrigCIDAC() getTrigMeV() getTrigRngMeV() --//
  for (XtalDiode xDiode; xDiode.isValid(); xDiode++) {
    const DiodeIdx diodeIdx(xtalIdx, xDiode);
    float tholdCIDAC;
    if (precalcCalibTool.getTrigCIDAC(diodeIdx, tholdCIDAC).isFailure())
      return StatusCode::FAILURE;

    //-- getTrigRngADC() --//
    RngNum rng;
    float tholdADC;
    if (precalcCalibTool.getTrigRngADC(diodeIdx, rng, tholdADC).isFailure())
      return StatusCode::SUCCESS;

    float testCIDAC;
    const RngIdx rngIdx(xtalIdx, xDiode.getFace(), rng);
    if (calCalibSvc.evalCIDAC(rngIdx, tholdADC, testCIDAC).isFailure())
      return StatusCode::FAILURE;
    
    if (!smart_compare(tholdCIDAC, testCIDAC, MAX_SPLINE_DIFF)) {
      MsgStream msglog(m_msgSvc, "test_PrecalcCalibTool");   
      msglog << MSG::ERROR << "Invalid trigger threshold: " << diodeIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    //-- getTrigMeV() test
    float trigMeV;
    if (precalcCalibTool.getTrigMeV(diodeIdx, trigMeV).isFailure())
      return StatusCode::FAILURE;

    float faceSignal;
    if (calCalibSvc.evalFaceSignal(rngIdx, tholdADC, faceSignal).isFailure())
      return StatusCode::FAILURE;

    if (!smart_compare(trigMeV, faceSignal, MAX_SPLINE_DIFF)) {
      MsgStream msglog(m_msgSvc, "test_PrecalcCalibTool");   
      msglog << MSG::ERROR << "Invalid trigger threshold: " << diodeIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }
  }

  //-- getLacCIDAC() test --//
  for (FaceNum face; face.isValid(); face++) {
    const FaceIdx faceIdx(xtalIdx, face);
    float lacCIDAC;
    if (precalcCalibTool.getLacCIDAC(faceIdx, lacCIDAC).isFailure())
      return StatusCode::FAILURE;

    CalibData::CalTholdCI const*const tholdCICalib = calCalibSvc.getTholdCI(faceIdx);
    if (!tholdCICalib) {
      MsgStream msglog(m_msgSvc, "test_CalCalibSvc");   
      msglog << MSG::ERROR << "MISSING tholdCI: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    const float lacADC = tholdCICalib->getLAC()->getVal();

    float testCIDAC;
    if (calCalibSvc.evalCIDAC(RngIdx(xtalIdx,face,LEX8), lacADC, testCIDAC).isFailure())
      return StatusCode::FAILURE;

    if (!smart_compare(lacCIDAC, testCIDAC, MAX_SPLINE_DIFF)) {
      MsgStream msglog(m_msgSvc, "test_PrecalcCalibTool");   
      msglog << MSG::ERROR << "Invalid lac threshold: " << faceIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }
  }

  return StatusCode::SUCCESS;
}


StatusCode test_PrecalcCalibTool::testMissingXtal(const XtalIdx xtalIdx,
                                                  IPrecalcCalibTool &precalcCalibTool) {
  float tmp;

  /// all calls should return failure
  for (XtalRng xRng; xRng.isValid(); xRng++) {
    const RngIdx rngIdx(xtalIdx, xRng);
    if (precalcCalibTool.getPedSigCIDAC(rngIdx, tmp).isSuccess())
      return StatusCode::FAILURE;
  }

  for (XtalDiode xDiode; xDiode.isValid(); xDiode++) {
    const DiodeIdx diodeIdx(xtalIdx, xDiode);
    if (precalcCalibTool.getTrigCIDAC(diodeIdx, tmp).isSuccess())
      return StatusCode::FAILURE;

    if (precalcCalibTool.getTrigMeV(diodeIdx, tmp).isSuccess())
      return StatusCode::FAILURE;

    RngNum rng;
    if (precalcCalibTool.getTrigRngADC(diodeIdx, rng, tmp).isSuccess())
      return StatusCode::FAILURE;

    
  }

  for (FaceNum face; face.isValid(); face++) {
    const FaceIdx faceIdx(xtalIdx, face);
    if (precalcCalibTool.getLacCIDAC(faceIdx, tmp).isSuccess())
      return StatusCode::FAILURE;
  }
  
  return StatusCode::SUCCESS;
}


StatusCode test_PrecalcCalibTool::verify(IPrecalcCalibTool &precalcCalibTool,
                                         ICalCalibSvc &calCalibSvc,
                                         const CalXtalResponse::TestCfg &testCfg,
                                         const TwrSet &twrSet) {
  MsgStream msglog(m_msgSvc, "test_PrecalcCalibTool");   

  // xtal loop
  for (TestCfg::XtalList::const_iterator xtalIt(testCfg.testXtals.begin());
       xtalIt != testCfg.testXtals.end();
       xtalIt++) {

    if (twrSet.find(xtalIt->getTwr()) == twrSet.end()) {
      if (testMissingXtal(*xtalIt,
                          precalcCalibTool).isFailure())
        return StatusCode::FAILURE;
    } else
      if (testXtal(*xtalIt,
                   calCalibSvc,
                   precalcCalibTool).isFailure())
        return StatusCode::FAILURE;
  } // xtal loop

  return StatusCode::SUCCESS;
}
                                                   


