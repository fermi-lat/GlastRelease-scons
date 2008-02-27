// $Header$
/** @file 
    @author Z.Fewtrell


*/

// LOCAL INCLUDES
#include "test_CalDiagnosticTool.h"
#include "CalXtalResponse/ICalDiagnosticTool.h"
#include "../CalCalib/IPrecalcCalibTool.h"
#include "CalXtalResponse/ICalSignalTool.h"
#include "CalXtalResponse/ICalTrigTool.h"
#include "test_util.h"
#include "TestCfg.h"

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"

// EXTLIB INCLUDES
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"

// STD INCLUDES
#include <memory>
#include <set>

using namespace CalUtil;
using namespace std;
using namespace CalXtalResponse;

/// loop thruogh & check cal diagnostic word for each tower & layer
StatusCode test_CalDiagnosticTool::verify(ICalSignalTool &calSignalTool,
                                          ICalTrigTool &calTrigTool,
                                          IPrecalcCalibTool &precalcCalibTool,
                                          ICalDiagnosticTool &calDiagnosticTool,
                                          const CalXtalResponse::TwrSet &twrSet) {

  for (TwrNum twr; twr.isValid(); twr++)
    for (LyrNum lyr; lyr.isValid(); lyr++)
      // CASE 1: active TEM
      if (twrSet.find(twr) != twrSet.end()) {
        if (testSingleDiag(calSignalTool,
                           calTrigTool,
                           precalcCalibTool,
                           calDiagnosticTool,
                           twr,
                           lyr).isFailure())
          return StatusCode::FAILURE;
      }
      // CASE 2: uninstalled TEM
      else {
        if (testEmptyDiag(calDiagnosticTool,
                          twr,
                          lyr).isFailure())
          return StatusCode::FAILURE;
      }

  return StatusCode::SUCCESS;
}

StatusCode test_CalDiagnosticTool::testSingleDiag(ICalSignalTool &calSignalTool,
                                                  ICalTrigTool &calTrigTool,
                                                  IPrecalcCalibTool &precalcCalibTool,
                                                  ICalDiagnosticTool &calDiagnosticTool,
                                                  const TwrNum twr,
                                                  const LyrNum lyr) {

  /// retrieve diagnostic bits
  auto_ptr<LdfEvent::CalDiagnosticData> diagData(calDiagnosticTool.getDiagnosticData(twr,lyr));
  if (diagData.get() == 0)
    return StatusCode::FAILURE;

  const unsigned dataWord = diagData->dataWord();

  /// check LAC bits
  for (FaceNum face; face.isValid(); face++)
    for (ColNum col; col.isValid(); col++) {
      float signalLevel;
      if (calSignalTool.getDiodeSignal(DiodeIdx(twr,lyr,col,face,LRG_DIODE),
                                       signalLevel).isFailure())
        return StatusCode::FAILURE;

      float thold;
      if (precalcCalibTool.getLacCIDAC(FaceIdx(twr,lyr,col,face),
                                       thold).isFailure())
        return StatusCode::FAILURE;
      
      /// expected lac value
      const bool expectedLAC = signalLevel > thold;

      // 1st bit in word for current GCRC
      const unsigned short faceFirstBit = (face == POS_FACE) ? 0 : 16;

      /// lac bit out of diagnostic word
      const bool diagLAC = dataWord & 1 << (faceFirstBit + col.val());

      if (diagLAC != expectedLAC) {
        MsgStream msglog(m_msgSvc, "test_CalDiagnosticTool"); 
        msglog << MSG::ERROR << "invalid diagnostic LAC bit " << FaceIdx(twr,lyr,col,face).toStr() << endreq;
        return StatusCode::FAILURE;
      }
    }

  /// check Trigger bits (check against 'or' of crystals all on one face
  for (FaceNum face; face.isValid(); face++)
    for (DiodeNum diode; diode.isValid(); diode++)  {
      bool expectedBit = false;

      for (ColNum col; col.isValid(); col++) {
        bool testBit;
        if (calTrigTool.getTriggerBit(DiodeIdx(twr,lyr,col,face,diode), testBit).isFailure())
          return StatusCode::FAILURE;
        
        expectedBit |= testBit;

        if (testBit) // no need to check any more columns if we found
                     // a raised trigger bit
          break;
      }

      // 1st bit in word for current GCRC
      const unsigned short faceFirstBit = (face == POS_FACE) ? 0 : 16;

      const bool diagBit = dataWord & 1 << (faceFirstBit + ColNum::N_VALS + diode.val());

      if (diagBit != expectedBit) {
        MsgStream msglog(m_msgSvc, "test_CalDiagnosticTool"); 
        msglog << MSG::ERROR << "invalid diagnostic trigger bit " << twr.toStr() 
               << " " << lyr.toStr() 
               << " " << face.toStr()
               << " " << diode.toStr() 
               << endreq;
        return StatusCode::FAILURE;
      }
      
      
    }

  return StatusCode::SUCCESS;
}

StatusCode test_CalDiagnosticTool::testEmptyDiag(ICalDiagnosticTool &calDiagnosticTool,
                                                 const CalUtil::TwrNum twr,
                                                 const CalUtil::LyrNum lyr) {
  
  /// retrieve diagnostic bits
  auto_ptr<LdfEvent::CalDiagnosticData> diagData(calDiagnosticTool.getDiagnosticData(twr,lyr));
  if (diagData.get() != 0)
    return StatusCode::FAILURE;

  return StatusCode::SUCCESS;
}

StatusCode test_CalDiagnosticTool::verifyTDS(LdfEvent::DiagnosticData *diagTds,
                                             ICalDiagnosticTool &calDiagnosticTool,
                                             const TestCfg &testCfg,
                                             const CalXtalResponse::TwrSet &twrSet) {
  /// check that we got diag data i.f.f. we asked for it.
  if ((diagTds !=0) != testCfg.createDiagnosticData) 
    return StatusCode::FAILURE;

  // quit if diag not enabled
  if (!testCfg.createDiagnosticData)
    return StatusCode::SUCCESS;

  // check number of diagnostic entries (8 per enabled TEM)
  const unsigned numCalDiag = diagTds->getNumCalDiagnostic();
  const unsigned numDiagExpected = twrSet.size()*LyrNum::N_VALS;
  if (numCalDiag != numDiagExpected) {
    MsgStream msglog(m_msgSvc, "test_CalDiagnosticTool"); 
    msglog << MSG::ERROR << "wrong # of calDiag "
           << numCalDiag << " in TDS (expected "
           << numDiagExpected << ")" << endreq;
    return StatusCode::FAILURE;
  }

  // keep track of all cal diagnostics retrieved so far
  CalVec<LyrIdx, bool> lyrSet;

  for (unsigned short idx = 0;
       idx < numCalDiag;
       idx++) {
    LdfEvent::CalDiagnosticData calDiagTds = diagTds->getCalDiagnosticByIndex(idx);

    TwrNum twr(calDiagTds.tower());
    LyrNum lyr(calDiagTds.layer());
    LyrIdx lyrIdx(twr,lyr);

    /// check for duplicate twr,lyr ident
    if (lyrSet[lyrIdx]) {
      MsgStream msglog(m_msgSvc, "test_CalDiagnosticTool"); 
      msglog << MSG::ERROR << "Duplicate CalDiag in TDS: " << twr.toStr() << " " << lyr.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    lyrSet[lyrIdx] = true;
    
    /// compare against answer from (already verified) cal diagnostic tool
    const unsigned diagWordTDS = calDiagTds.dataWord();
    
    /// retrieve diagnostic bits
    auto_ptr<LdfEvent::CalDiagnosticData> diagData(calDiagnosticTool.getDiagnosticData(twr,lyr));
    if (diagData.get() == 0)
      return StatusCode::FAILURE;

    const unsigned diagWordTool = diagData->dataWord();
    if (diagWordTDS != diagWordTool) {
      MsgStream msglog(m_msgSvc, "test_CalDiagnosticTool"); 
      msglog << MSG::ERROR << "Mismatch cal diag word: " << diagWordTool << " vs TDS:" << diagWordTDS << endreq;
      return StatusCode::FAILURE;
    }


  }

  return StatusCode::SUCCESS;
}
