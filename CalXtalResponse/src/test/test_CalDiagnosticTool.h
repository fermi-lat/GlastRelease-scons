#ifndef test_CalDiagnosticTool_h
#define test_CalDiagnosticTool_h
//     $Header$

/** @file
    @author Z.Fewtrell
    unit test for CalDiagnosticTool class
*/

// LOCAL INCLUDES
#include "TestCfg.h"

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"

// EXTLIB INCLUDES
#include "GaudiKernel/StatusCode.h"

// STD INCLUDES


class ICalSignalTool;
class IPrecalcCalibTool;
class ICalTrigTool;
class ICalDiagnosticTool;
class IMessageSvc;

namespace LdfEvent {
  class DiagnosticData;
};

/** unit test for CalDiagnosticTool
 */
class test_CalDiagnosticTool {
public:

  /// Gaudi initialize() phase intializations
  StatusCode initialize(IMessageSvc *msgSvc) {
    m_msgSvc = msgSvc;

    return StatusCode::SUCCESS;
  }

  /// verify given CalDiagnosticTool object w/ given valid TEM list & test calibration set
  StatusCode verify(ICalSignalTool &calSignalTool,
                    ICalTrigTool &calTrigTool,
                    IPrecalcCalibTool &precalcCalibTool,
                    ICalDiagnosticTool &calDiagnosticTool,
                    const CalXtalResponse::TwrSet &twrSet);

  /// verify that Cal Diagnostic data in TDS is correct
  StatusCode verifyTDS(LdfEvent::DiagnosticData *diagTds,
                       ICalDiagnosticTool &calDiagnosticTool,
                       const CalXtalResponse::TestCfg &testCfg,
                       const CalXtalResponse::TwrSet &twrSet);

private:
  /// test all bits for single cal diagnostic word 
  StatusCode testSingleDiag(ICalSignalTool &calSignalTool,
                            ICalTrigTool &calTrigTool,
                            IPrecalcCalibTool &precalcCalibTool,
                            ICalDiagnosticTool &calDiagnosticTool,
                            const CalUtil::TwrNum twr,
                            const CalUtil::LyrNum lyr);

  /// test correct response for uninstalled tem.
  StatusCode testEmptyDiag(ICalDiagnosticTool &calDiagnosticTool,
                           const CalUtil::TwrNum twr,
                           const CalUtil::LyrNum lyr);

  
  /// Gleam logging service
  IMessageSvc *m_msgSvc;
};

#endif // test_CalDiagnosticTool_h
