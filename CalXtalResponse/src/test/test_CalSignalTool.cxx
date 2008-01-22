// $Header$
/** @file 
    @author Z.Fewtrell


*/

// LOCAL INCLUDES
#include "test_CalSignalTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"
#include "CalXtalResponse/ICalSignalTool.h"
#include "test_util.h"
#include "TestCfg.h"

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "CalibData/Cal/CalMevPerDac.h"
#include "CalUtil/CalArray.h"

// EXTLIB INCLUDES

// STD INCLUDES
#include <set>
#include <iterator>
#include <cmath>

using namespace CalUtil;
using namespace std;
using namespace CalXtalResponse;



/// verify given CalSignalTool object w/ given valid TEM list & test calibration set
StatusCode test_CalSignalTool::verify(ICalSignalTool &calSignalTool,
                                      ICalCalibSvc &calCalibSvc,
                                      const CalXtalResponse::TestCfg &testCfg,
                                      const CalXtalResponse::TwrSet &twrSet) {
  if (testCalRelationMap(calSignalTool,
                         testCfg,
                         twrSet).isFailure())
    return StatusCode::FAILURE;
  
  typedef set<XtalIdx> XtalSet;
  XtalSet xtalSet(testCfg.testXtals.begin(), testCfg.testXtals.end());
  
  for (XtalIdx xtalIdx; xtalIdx.isValid(); xtalIdx++) {
    // any crystal without hit should be empty as we have no noise enabled
    if (xtalSet.find(xtalIdx) == xtalSet.end() || twrSet.find(xtalIdx.getTwr()) == twrSet.end()) {
      if (checkEmptyXtal(xtalIdx, calSignalTool).isFailure())
        return StatusCode::FAILURE;
    }
    else if (testXtal(xtalIdx, calSignalTool, calCalibSvc, testCfg).isFailure())
      return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}


/** test all signal level for all diodes on given crystal
    Algorithm:
    - calculate mean CIDAC & CIDAC asymmetry from CalSignalTool
    - using mev & asymmetry calculations, determine energy deposit
    position and magnitude
    - compare against test case configuration in testCfg
*/
StatusCode test_CalSignalTool::testXtal(const CalUtil::XtalIdx xtalIdx,
                                        ICalSignalTool &calSignalTool,
                                        ICalCalibSvc &calCalibSvc,
                                        const CalXtalResponse::TestCfg &testCfg) {

  /// retrieve signal for all diodes in xtal
  CalArray<XtalDiode, float> cidac(0);
  for (XtalDiode xDiode;
       xDiode.isValid();
       xDiode++) {
    if (calSignalTool.getDiodeSignal(DiodeIdx(xtalIdx, xDiode), cidac[xDiode]).isFailure()) {
      MsgStream msglog(m_msgSvc, "test_CalSignalTool");   
      msglog << MSG::ERROR << "getDiodeSignal failure: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }
  }

  /// retrieve calibration constants
  CalibData::CalMevPerDac const*const mpdCalib = calCalibSvc.getMPD(xtalIdx);
  if (!mpdCalib) {
    MsgStream msglog(m_msgSvc, "test_CalSignalTool");   
    msglog << MSG::ERROR << "calib failure: " << xtalIdx.toStr() << endreq;
    return StatusCode::FAILURE;
  }

  CalArray<DiodeNum, float> mpd(0);
  mpd[LRG_DIODE] = mpdCalib->getBig()->getVal();
  mpd[SM_DIODE] = mpdCalib->getSmall()->getVal();

  /// check energy level & asymmetry for both diodes
  for (DiodeNum diode; diode.isValid(); diode++) {
    const float meanDAC = sqrt(cidac[XtalDiode(POS_FACE, diode)] *
                               cidac[XtalDiode(NEG_FACE, diode)]);

    const float meV = meanDAC*mpd[diode];
    const float testMeV = testCfg.totalEnergy();

    float pctTolerance = MAX_RECON_DIFF;
    // high curvature in asymmetry model introduces error when summing hits
    if (testCfg.xtalHits.size() > 1)
      pctTolerance = MAX_ASYM_DIFF;

    if (!smart_compare(meV, testMeV, pctTolerance)) {
      MsgStream msglog(m_msgSvc, "test_CalSignalTool");   
      msglog << MSG::ERROR << "bad signal level: " << xtalIdx.toStr() << " " << diode.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    const float asym = log(cidac[XtalDiode(POS_FACE, diode)] / 
                           cidac[XtalDiode(NEG_FACE, diode)]);
    
    float posMM = 0;
    if (calCalibSvc.evalPos(xtalIdx,
                            AsymType(diode,diode),
                            asym,
                            posMM).isFailure())
      return StatusCode::FAILURE;
    
    /// convert xtal (-167...167) position in mm to relative position (0...1)
    const float relPos = xtalMM2RelPos(posMM);

    if (abs_diff(relPos, testCfg.weightedPos()) > MAX_RECON_DIFF*csiMM) {
      MsgStream msglog(m_msgSvc, "test_CalSignalTool");   
      msglog << MSG::ERROR << "bad xtal asym: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }
  }
  
                                        
  return StatusCode::SUCCESS;
}

/** Algorithm:
    - check correct number of xtals registered (account for any
    missing tower modules)
    - check correct number of hits in each xtal against testCfg
    - check total ene in all MC hits for xtal match total ene in testCfg
 */
StatusCode test_CalSignalTool::testCalRelationMap(ICalSignalTool &calSignalTool,
                                                  const CalXtalResponse::TestCfg &testCfg,
                                                  const CalXtalResponse::TwrSet &twrSet) {

  const ICalSignalTool::CalRelationMap *calRelMap = calSignalTool.getCalRelationMap();
  if (!calRelMap) {
    MsgStream msglog(m_msgSvc, "test_CalSignalTool");   
    msglog << MSG::ERROR << "missing cal relation table: " << endreq;
    return StatusCode::FAILURE;
  }

  /// check that right number of hits have been registered
  /// find how many xtals are valid (i.e. the towers are installed)
  const unsigned short nValidXtals = testCfg.nValidXtals(twrSet);
  if (calRelMap->size() != testCfg.xtalHits.size()*nValidXtals) {
    MsgStream msglog(m_msgSvc, "test_CalSignalTool");   
    msglog << MSG::ERROR << "wrong # of cal relations" << endreq;
    return StatusCode::FAILURE;
  }

  /// now check for individual hits
  for (TestCfg::XtalList::const_iterator xtalIt(testCfg.testXtals.begin());
       xtalIt != testCfg.testXtals.end();
       xtalIt++) {
    
    /// skip xtals in disabled towers
    if (twrSet.find(xtalIt->getTwr()) == twrSet.end())
      continue;

    typedef ICalSignalTool::CalRelationMap::const_iterator CalRelationIt;
    pair<CalRelationIt, CalRelationIt> xtalHitMatches(calRelMap->equal_range(*xtalIt));

    /// check right # of hits
    if (distance(xtalHitMatches.first, xtalHitMatches.second) != static_cast<int>(testCfg.xtalHits.size())) {
      MsgStream msglog(m_msgSvc, "test_CalSignalTool");   
      msglog << MSG::ERROR << "wrong # of xtal hits: " << xtalIt->toStr() << endreq;
      return StatusCode::FAILURE;
    }

    /// check correct total energy & xtal match
    float totalMeV = 0;
    for (CalRelationIt it(xtalHitMatches.first);
         it != xtalHitMatches.second;
         it++) {
      const XtalIdx hitIdx(idents::CalXtalId(it->second->volumeID()));
      if (*xtalIt != hitIdx) {
        MsgStream msglog(m_msgSvc, "test_CalSignalTool");   
        msglog << MSG::ERROR << "CalRelTable xtalIdx mismatch" << endreq;
        return StatusCode::FAILURE;
      }

      totalMeV += it->second->totalEnergy();
    }

    if (!smart_compare(totalMeV,testCfg.totalEnergy(), MAX_RECON_DIFF)) {
      MsgStream msglog(m_msgSvc, "test_CalSignalTool");   
      msglog << MSG::ERROR << "CalRelTable bad total energy: " << xtalIt->toStr() << endreq;
      return StatusCode::FAILURE;
    }
  }
    
  return StatusCode::SUCCESS;
}


StatusCode test_CalSignalTool::checkEmptyXtal(const CalUtil::XtalIdx xtalIdx,
                                              ICalSignalTool &calSignalTool) {
  for (XtalDiode xDiode;
       xDiode.isValid();
       xDiode++) {
    float tmp = -1;
    if (calSignalTool.getDiodeSignal(DiodeIdx(xtalIdx, xDiode), tmp).isFailure()) {
      MsgStream msglog(m_msgSvc, "test_CalSignalTool");   
      msglog << MSG::ERROR << "getDiodeSignal failure: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }

    if (tmp != 0) {
      MsgStream msglog(m_msgSvc, "test_CalSignalTool");   
      msglog << MSG::ERROR << "'empty' xtal failure: " << xtalIdx.toStr() << endreq;
      return StatusCode::FAILURE;
    }
  }

  return StatusCode::SUCCESS;
}
