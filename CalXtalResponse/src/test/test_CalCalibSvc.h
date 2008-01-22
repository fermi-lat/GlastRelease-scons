#ifndef test_CalCalibSvc_h
#define test_CalCalibSvc_h
//     $Header$

/** @file
    @author Z.Fewtrell
    
   unit test for CalCalibSvc class
*/

// LOCAL INCLUDES
#include "TestCfg.h"

// GLAST INCLUDES
#include "CalUtil/SimpleCalCalib/CalPed.h"
#include "CalUtil/SimpleCalCalib/CalAsym.h"
#include "CalUtil/SimpleCalCalib/CalMPD.h"
#include "CalUtil/SimpleCalCalib/CIDAC2ADC.h"
#include "CalUtil/SimpleCalCalib/TholdCI.h"

// EXTLIB INCLUDES
#include "GaudiKernel/StatusCode.h"

// STD INCLUDES
#include <set>

class IMessageSvc;
namespace CalXtalResponse {
  class TestCfg;
}
class ICalCalibSvc;

/** unit test for CalCalibSvc class
 */
class test_CalCalibSvc {
public:
  /** complete set of Cal Calibration constants, stored in arrays &
      loaded from simple txt files.  Good for testing ouptut of
      'official' calibration data stream.
   */
  class TestCalibSet {
  public:
    TestCalibSet() {};
    
    TestCalibSet(const std::string &pedTXTPath,
                 const std::string &cidac2adcTXTPath,
                 const std::string &asymTXTPath,
                 const std::string &mpdTXTPath,
                 const std::string &tholdCITXTPath);

    
    
    CalUtil::CalPed m_calPed;
    CalUtil::CIDAC2ADC m_cidac2adc;
    CalUtil::CalAsym m_calAsym;
    CalUtil::CalMPD m_calMPD;
    CalUtil::TholdCI m_calTholdCI;
  };

  /// Gaudi initialize() phase intializations
  StatusCode initialize(IMessageSvc *msgSvc) {
    m_msgSvc = msgSvc;

    return StatusCode::SUCCESS;
  }

  
  /// verify given CalCalibSvc object w/ given valid TEM list & test calibration set
  StatusCode verify(ICalCalibSvc &calCalibSvc,
                    const test_CalCalibSvc::TestCalibSet &calibSet,
                    const CalXtalResponse::TestCfg &testCfg,
					const CalXtalResponse::TwrSet &twrSet);

  

private:
  
  /// test all calibrations on given xtal
  StatusCode testXtal(const CalUtil::XtalIdx xtalIdx,
                      ICalCalibSvc &calCalibSvc,
                      const TestCalibSet &calibSet,
                      const CalXtalResponse::TwrSet &twrSet);

  /// test that missing xtal (from uninstalled tower) behaves
  /// properly (i.e. always returns null)
  StatusCode testMissingXtal(const CalUtil::XtalIdx xtalIdx,
                             ICalCalibSvc &calCalibSvc);


  /// verify pedestal calibration for single xtal
  StatusCode testPed(const CalUtil::XtalIdx xtalIdx,
                     ICalCalibSvc &calCalibSvc,
                     const TestCalibSet &calibSet);

  /// verify intNonlin calibration for single xtal
  StatusCode testINL(const CalUtil::XtalIdx xtalIdx,
                     ICalCalibSvc &calCalibSvc,
                     const TestCalibSet &calibSet);

  /// verify asymmetry calibration for single xtal
  StatusCode testAsym(const CalUtil::XtalIdx xtalIdx,
                      ICalCalibSvc &calCalibSvc,
                      const TestCalibSet &calibSet);

  /// verify mevPerDAC calibration for single xtal
  StatusCode testMPD(const CalUtil::XtalIdx xtalIdx,
                     ICalCalibSvc &calCalibSvc,
                     const TestCalibSet &calibSet);

  /// verify threshold calibrations for single xtal
  StatusCode testTholdCI(const CalUtil::XtalIdx xtalIdx,
                         ICalCalibSvc &calCalibSvc,
                         const TestCalibSet &calibSet);

  /// verify faceSignal (adc->MeV) calculation for single crystal
  StatusCode testFaceSignal(const CalUtil::XtalIdx xtalIdx,
                            ICalCalibSvc &calCalibSvc);

  /// Gleam logging object
  IMessageSvc *m_msgSvc;
};

#endif // test_CalCalibSvc_h
