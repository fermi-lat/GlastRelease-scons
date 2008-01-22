#ifndef test_CalSignalTool_h
#define test_CalSignalTool_h
//     $Header$

/** @file
    @author Z.Fewtrell
    unit test for CalSignalTool class
*/

// LOCAL INCLUDES
#include "TestCfg.h"

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"

// EXTLIB INCLUDES
#include "GaudiKernel/StatusCode.h"

// STD INCLUDES
#include <vector>

class IMessageSvc;
namespace CalXtalResponse {
  class TestCfg;
}
class ICalSignalTool;
class ICalCalibSvc;

/** unit test for CalSignalTool
 */
class test_CalSignalTool {
public:

  /// Gaudi initialize() phase intializations
  StatusCode initialize(IMessageSvc *msgSvc) {
    m_msgSvc = msgSvc;

    return StatusCode::SUCCESS;
  }

  /// verify given CalSignalTool object w/ given valid TEM list & test calibration set
  StatusCode verify(ICalSignalTool &calSignalTool,
                    ICalCalibSvc &calCalibSvc,
                    const CalXtalResponse::TestCfg &testCfg,
                    const CalXtalResponse::TwrSet &twrSet);

private:
  
  /// test all calibrations on given xtal
  StatusCode testXtal(const CalUtil::XtalIdx xtalIdx,
                      ICalSignalTool &calSignalTool,
                      ICalCalibSvc &calCalibSvc,
                      const CalXtalResponse::TestCfg &testCfg);

  /// check xtalId->Mc hit relation map
  StatusCode testCalRelationMap(ICalSignalTool &calSignalTool,
                                const CalXtalResponse::TestCfg &testCfg,
                                const CalXtalResponse::TwrSet &twrSet);

  /// check that response for uninstalled crystal (i.e. missing tower
  /// module) is correct
  StatusCode checkEmptyXtal(const CalUtil::XtalIdx xtalIdx,
                            ICalSignalTool &calSignalTool);

  /// Gleam logging service
  IMessageSvc *m_msgSvc;
};

#endif // test_CalSignalTool_h
