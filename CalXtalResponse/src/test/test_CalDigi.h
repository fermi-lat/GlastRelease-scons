#ifndef test_CalDigi_h
#define test_CalDigi_h
//     $Header$

/** @file
    @author Z.Fewtrell
     unit test for CalDigi class
*/

// LOCAL INCLUDES
#include "TestCfg.h"

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"
#include "Event/Digi/CalDigi.h"

// EXTLIB INCLUDES
#include "GaudiKernel/StatusCode.h"

// STD INCLUDES
#include <vector>

class IMessageSvc;
namespace CalXtalResponse {
  class TestCfg;
}
class ICalCalibSvc;
class ICalSignalTool;


/** \brief unit test for CalDigiAlg results 
 */
class test_CalDigi {
public:

  /// Gaudi initialize() phase intializations
  StatusCode initialize(IMessageSvc *msgSvc) {
    m_msgSvc = msgSvc;

    return StatusCode::SUCCESS;
  }

  /// verify given CalDigi object w/ given valid TEM list & test calibration set
  StatusCode verify(ICalSignalTool &calSignalTool,
                    ICalCalibSvc &calCalibSvc,
                    const CalXtalResponse::TestCfg &testCfg,
                    const CalXtalResponse::TwrSet &twrSet,
                    const Event::CalDigiCol &calDigiCol);
  
private:

  /// verify individual cal digi
  StatusCode verifyXtal(ICalSignalTool &calSignalTool,
                        ICalCalibSvc &calCalibSvc,
                        const CalXtalResponse::TestCfg &testCfg,
                        const Event::CalDigi &calDigi);

  /// check for duplicated CalXtalId in digi collection
  StatusCode checkDuplicateXtals(const Event::CalDigiCol &calDigiCol);

  /// check number of readouts in CalDigi object against current
  /// software configuration
  StatusCode checkNReadouts(const CalXtalResponse::TestCfg &testCfg,
                            const Event::CalDigi &calDigi);

  /// check best range selection in CalDigi
  StatusCode checkBestRange(ICalSignalTool &calSignalTool,
                            ICalCalibSvc &calCalibSvc,
                            const Event::CalDigi &calDigi
                            );

  /// check readout range order in CalDigi object.
  StatusCode checkRangeOrder(const Event::CalDigi &calDigi);

  /// check ADC readout values in CalDigi
  StatusCode checkADC(ICalSignalTool &calSignalTool,
                      ICalCalibSvc &calCalibSvc,
                      const Event::CalDigi &calDigi);

  StatusCode checkMissingXtals(const Event::CalDigiCol &calDigiCol,
                               const CalXtalResponse::TestCfg &testCfg);


  /// Gleam loggin service
  IMessageSvc *m_msgSvc;
};

#endif // test_CalSignalTool_h
