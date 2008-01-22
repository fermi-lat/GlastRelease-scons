#ifndef test_CalTrigTool_h
#define test_CalTrigTool_h
//     $Header$

/** @file
    @author Z.Fewtrell
   unit test for CalTrigTool class
*/

// LOCAL INCLUDES
#include "TestCfg.h"

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"

// EXTLIB INCLUDES
#include "GaudiKernel/StatusCode.h"

// STD INCLUDES

class IMessageSvc;
namespace CalXtalResponse {
  class TestCfg;
}
class ICalTrigTool;
class ICalCalibSvc;
class ICalSignalTool;
class IPrecalcCalibTool;

/** unit test for CalTrigTool
 */
class test_CalTrigTool {
public:

  /// Gaudi initialize() phase intializations
  StatusCode initialize(IMessageSvc *msgSvc) {
    m_msgSvc = msgSvc;

    return StatusCode::SUCCESS;
  }

  /// verify given CalTrigTool object w/ given valid TEM list & test calibration set
  StatusCode verify(ICalSignalTool &calSignalTool,
                    const CalXtalResponse::TwrSet &twrSet,
                    ICalTrigTool &calTrigTool,
                    IPrecalcCalibTool &precalcCalibTool,
					const float tolerance);
  
private:

  /// verify individual cal digi
  StatusCode verifyDiode(ICalSignalTool &calSignalTool,
                         ICalTrigTool &calTrigTool,
                         const CalUtil::DiodeIdx diodeIdx,
                         IPrecalcCalibTool &precalcCalibTool,
					     const float tolerance);

  /// verify per-tower trigger bits
  StatusCode verifyTriggerVecs(ICalSignalTool &calSignalTool,
                               ICalTrigTool &calTrigTool,
                               const CalXtalResponse::TwrSet &twrSet,
								IPrecalcCalibTool &precalcCalibTool,
							   const float tolerance);

  /// verify individual cal digi
  StatusCode verifyEmptyDiode(ICalTrigTool &calTrigTool,
                              const CalUtil::DiodeIdx diodeIdx);

  /// verify tower wide trigger bit
  StatusCode verifyTwrTrigger(const CalUtil::TwrNum twr,
                              const CalUtil::DiodeNum diode,
                              ICalSignalTool &calSignalTool,
                              IPrecalcCalibTool &precalcCalibTool,
                              const bool trigBit,
                              float tolerancePct);

  /// verify 16 bit (1 bit per tower) trigger vectors 
  StatusCode verifyTriggerVecs(ICalSignalTool &calSignalTool,
                               ICalTrigTool &calTrigTool,
                               const CalXtalResponse::TwrSet &twrSet,
                               IPrecalcCalibTool &precalcCalibTool);


  /// for logging
  IMessageSvc *m_msgSvc;
};

#endif // test_CalSignalTool_h
