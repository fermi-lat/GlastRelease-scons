#ifndef test_CalXtalRec_h
#define test_CalXtalRec_h
//     $Header$

/** @file
    @author Z.Fewtrell
    unit test for CalXtalRec class
*/

// LOCAL INCLUDES
#include "TestCfg.h"

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"
#include "Event/Digi/CalDigi.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"

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
class IPrecalcCalibTool;
class IGlastDetSvc;


/** unit test for CalXtalRecAlg results 
 */
class test_CalXtalRecAlg {
 public:

  /// Gaudi initialize() phase intializations
  StatusCode initialize(IMessageSvc *msgSvc,
                        IGlastDetSvc *detSvc);

  /// verify given CalXtalRec object w/ given valid TEM list &/** unit test for calibration set
  StatusCode verify(ICalSignalTool &calSignalTool,
                    ICalCalibSvc &calCalibSvc,
		    IPrecalcCalibTool &precalcCalibTool,
                    const CalXtalResponse::TestCfg &testCfg,
                    const CalXtalResponse::TwrSet &twrSet,
		    const Event::CalDigiCol &calDigiCol,
                    const Event::CalXtalRecCol &calXtalRecCol);
  
 private:
  /// verify individual cal digi
  StatusCode verifyXtal(ICalCalibSvc &calCalibSvc,
                        const CalXtalResponse::TestCfg &testCfg,
			const Event::CalDigi &calDigi,
                        const Event::CalXtalRecData &calXtalRecData);

  /// compare recon readout vs corresponding digi readout
  StatusCode verifyRangeReadout(const CalUtil::XtalIdx xtalIdx,
                                ICalCalibSvc &calCalibSvc,
                                const Event::CalDigi::CalXtalReadout &digiRO,
				const Event::CalXtalRecData::CalRangeRecData &reconRO);

  /// compare recon readout vs MC truth.
  StatusCode verifyRangeRecon(const CalUtil::XtalIdx xtalIdx,
                              ICalCalibSvc &calCalibSvc,
                              const CalXtalResponse::TestCfg &testCfg,
							  const Event::CalDigi::CalXtalReadout &digiRO,
							  const Event::CalXtalRecData::CalRangeRecData &reconRO);

  /** \brief convert longitudinal xtalPos along cal xtal to 3d point
                                  in LAT geometry space

    \param pXtal ouput position vector
    \param xtalPos input longitudinal position in mm from center of
    xtal
  */
  Point pos2Point(const CalUtil::XtalIdx xtalIdx, 
                  const float xtalPosMM) const;

  /// Gleam logging service
  IMessageSvc *m_msgSvc;

  /// the value of fLATObjects field, defining LAT towers 
  int m_eLATTowers; 
  
  /// the value of fTowerObject field, defining calorimeter module 
  int m_eTowerCAL;
  /// the value of fCellCmp field defining CsI crystal
  int m_eXtal;      
  /// number of geometric segments per Xtal
  int m_nCsISeg; 

  /// len of one Cal xtal
  float m_csiLength;

  /// used for constants & conversion routines.
  IGlastDetSvc* m_detSvc;
};

#endif // test_CalSignalTool_h

