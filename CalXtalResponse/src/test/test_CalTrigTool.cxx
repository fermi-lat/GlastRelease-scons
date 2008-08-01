// $Header$
/** @file
    @author Z.Fewtrell
*/

// LOCAL INCLUDES
#include "test_CalTrigTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"
#include "CalXtalResponse/ICalSignalTool.h"
#include "CalXtalResponse/ICalTrigTool.h"
#include "test_util.h"
#include "TestCfg.h"
#include "../CalCalib/IPrecalcCalibTool.h"

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"


// EXTLIB INCLUDES


// STD INCLUDES


using namespace CalUtil;
using namespace std;
using namespace CalXtalResponse;

/**\param tolerance percent tolerance allowed for trigger evaluation
   
   Algorithm:
   - verify trigger bit for each diode in Cal
   - verify CalTriggerVecs (per tower trigger bits)
 */
StatusCode test_CalTrigTool::verify(ICalSignalTool &calSignalTool,
                                    const CalXtalResponse::TwrSet &twrSet,
                                    ICalTrigTool &calTrigTool,
                                    IPrecalcCalibTool &precalcCalibTool,
                                    const float tolerance) {


  for (DiodeIdx diodeIdx; diodeIdx.isValid(); diodeIdx++) {
    /// CASE 0: active diode
    if (twrSet.find(diodeIdx.getTwr()) != twrSet.end()) {
      if (verifyDiode(calSignalTool,
                      calTrigTool,
                      diodeIdx,
                      precalcCalibTool,
                      tolerance).isFailure())
        return StatusCode::FAILURE;
    }
    else // CASE 1: inactive diode
      if (verifyEmptyDiode(calTrigTool,
                           diodeIdx).isFailure())
        return StatusCode::FAILURE;
                     
  }

  if (verifyTriggerVecs(calSignalTool,
                        calTrigTool,
                        twrSet,
                        precalcCalibTool,
                        tolerance).isFailure())
    return StatusCode::FAILURE;
        
  

  return StatusCode::SUCCESS;
}

/** Algorithm:
   - estimate trigger response with signal level (calSignalTool)& trigger threshold
   calibration (precalcCalibTool)
   - verifgy that calTrigTool response matches estimated response w/in
   tolerance percent of threshold.
 */
StatusCode test_CalTrigTool::verifyDiode(ICalSignalTool &calSignalTool,
                                         ICalTrigTool &calTrigTool,
                                         const DiodeIdx diodeIdx,
                                         IPrecalcCalibTool &precalcCalibTool,
                                         const float tolerance
                                         ) {
  bool trigBit = false;
  if (calTrigTool.getTriggerBit(diodeIdx, trigBit).isFailure())
    return StatusCode::FAILURE;

  /// get trigger threshold in CIDAC
  float trigCIDAC = 0;
  if (precalcCalibTool.getTrigCIDAC(diodeIdx, trigCIDAC).isFailure())
    return StatusCode::FAILURE;

  
  /// get signal level
  float signal = 0;
  if (calSignalTool.getDiodeSignal(diodeIdx, signal).isFailure())
    return StatusCode::FAILURE;

  /// test trigger threshold with some margin of error.
  if (!trig_test_margin(signal, trigCIDAC, trigBit, tolerance*signal)) {
    MsgStream msglog(m_msgSvc, "test_CalTrigTool"); 
    msglog << MSG::ERROR << "invalid trigger bit " << diodeIdx.toStr() << endreq;
    return StatusCode::FAILURE;
  }


  /// still need to test caltrigger vector
  return StatusCode::SUCCESS;
}

StatusCode test_CalTrigTool::verifyEmptyDiode(ICalTrigTool &calTrigTool,
                                              const CalUtil::DiodeIdx diodeIdx) {
  bool trigBit = true;

  if (calTrigTool.getTriggerBit(diodeIdx, trigBit).isFailure())
    return StatusCode::FAILURE;

  /// trigger should always be false for inactive xtal
  if (trigBit != false)
    return StatusCode::FAILURE;

  
  return StatusCode::SUCCESS;
}


StatusCode test_CalTrigTool::verifyTriggerVecs(ICalSignalTool &calSignalTool,
                                               ICalTrigTool &calTrigTool,
                                               const CalXtalResponse::TwrSet &twrSet,
                                               IPrecalcCalibTool &precalcCalibTool,
                                               const float tolerance) {

  for (DiodeNum diode; diode.isValid(); diode++) {
    unsigned short triggerVec;

    if (calTrigTool.getCALTriggerVector(static_cast<idents::CalXtalId::DiodeType>(diode.val()), triggerVec).isFailure())
      return StatusCode::FAILURE;

    for (TwrNum twr; twr.isValid(); twr++) {
      const bool trigBit = triggerVec & (1 << twr.val());
      if (twrSet.find(twr) != twrSet.end()) { // DEFAULT CASE (tower bay full)
        if (verifyTwrTrigger(twr, diode, calSignalTool, precalcCalibTool, trigBit, tolerance).isFailure())
          return StatusCode::FAILURE;
      } else { // CASE 2: (tower bay empty)
        /// trigger bit better be false
        if (trigBit) {
          MsgStream msglog(m_msgSvc, "test_CalTrigTool"); 
          msglog << MSG::ERROR << "invalid trigger vector bit in empty tower bay" << twr.toStr() << endreq;
          return StatusCode::FAILURE;
        }
      }
    }
  }

  return StatusCode::SUCCESS;
}

/** \param twr selected tower
    \param diode selected diode size (SMALL or LARGE)
    \param trigBit trigger bit from calTrigTool

   loop through every crystal face in tower & find whether tower
   trigger is certainly true, certainly false, or unknown (within
   tolerance).

   compare against trigBit.

   if trigger is within tolerancePct, return success
 */
StatusCode test_CalTrigTool::verifyTwrTrigger(const CalUtil::TwrNum twr,
                                              const CalUtil::DiodeNum diode,
                                              ICalSignalTool &calSignalTool,
                                              IPrecalcCalibTool &precalcCalibTool,
                                              const bool trigBit,
                                              float tolerancePct) {
  typedef enum ThreeBit { ThreeFalse,
                          ThreeUnknown,
                          ThreeTrue } ;

  ThreeBit triggerStatus(ThreeFalse);


  /// loop through every crystal face in tower
  for (tFaceIdx faceIdx;
       faceIdx.isValid();
       faceIdx++) {
    const DiodeIdx diodeIdx(twr, 
                            faceIdx.getLyr(), 
                            faceIdx.getCol(), 
                            faceIdx.getFace(), 
                            diode);

    /// get diode signal level
    float signal;
    if (calSignalTool.getDiodeSignal(diodeIdx, signal).isFailure())
      return StatusCode::FAILURE;

    /// get trigger threshold
    float thold;
    if (precalcCalibTool.getTrigCIDAC(diodeIdx, thold).isFailure())
      return StatusCode::FAILURE;

    /// bit is certainly true.
    if (!trigBit && signal > thold + thold*tolerancePct) {
      triggerStatus = ThreeTrue;
      /// quick return 
      if (trigBit)
        return StatusCode::SUCCESS;
      else {
        MsgStream msglog(m_msgSvc, "test_CalTrigTool"); 
        msglog << MSG::ERROR << "missing trigger bit in CalTriggerVec" << diodeIdx.toStr() << endreq;
        return StatusCode::FAILURE;
      }
    } else if (signal > thold - thold*tolerancePct && triggerStatus == ThreeFalse) // CASE within margin of error
      triggerStatus = ThreeUnknown;
  }

  /// case where trigger bit could go either way
  if (triggerStatus == ThreeUnknown)
    return StatusCode::SUCCESS;

  /// better be false @ this point
  if (trigBit) {
    MsgStream msglog(m_msgSvc, "test_CalTrigTool"); 
    msglog << MSG::ERROR << "false positive trigger bit in CalTriggerVec" << twr.toStr() << endreq;
    return StatusCode::FAILURE;
  }
  
  return StatusCode::SUCCESS;
}
