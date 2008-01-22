#ifndef test_PrecalcCalibTool_h
#define test_PrecalcCalibTool_h
//     $Header$

/** @file
    @author Z.Fewtrell
   unit test for PrecalcCalibTool class
*/

// LOCAL INCLUDES
#include "TestCfg.h"

// GLAST INCLUDES

// EXTLIB INCLUDES
#include "GaudiKernel/StatusCode.h"

// STD INCLUDES

class IMessageSvc;
namespace CalXtalResponse {
  class TestCfg;
}
class IPrecalcCalibTool;

/** unit test for PrecalcCalibTool
 */
class test_PrecalcCalibTool {
public:
  /// Gaudi initialize() phase intializations
  StatusCode initialize(IMessageSvc *msgSvc) {
    m_msgSvc = msgSvc;

    return StatusCode::SUCCESS;
  }

  
  /// verify given PrecalcCalibTool object w/ given valid TEM list & test calibration set
  StatusCode verify(IPrecalcCalibTool &precalcCalibTool,
                                         ICalCalibSvc &calCalibSvc,
                                         const CalXtalResponse::TestCfg &testCfg,
										 const CalXtalResponse::TwrSet &twrSet);

  

private:
  
  /// test all calibrations on given xtal
  StatusCode testXtal(const CalUtil::XtalIdx xtalIdx,
                      ICalCalibSvc &calCalibSvc,
                      IPrecalcCalibTool &precalcCalibTool);

  /// test that missing xtal (from uninstalled tower) behaves
  /// properly (i.e. always returns null)
  StatusCode testMissingXtal(const CalUtil::XtalIdx xtalIdx,
                             IPrecalcCalibTool &precalcCalibTool);

  /// Gleam logging object
  IMessageSvc *m_msgSvc;
};

#endif // test_PrecalcCalibTool_h
