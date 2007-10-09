// $Header$

// Include files
// Gaudi system includes

// LOCAL
#include "CalXtalResponse/ICalCalibSvc.h"
#include "CalXtalResponse/IXtalRecTool.h"
#include "CalXtalResponse/IXtalDigiTool.h"
#include "CalXtalResponse/ICalSignalTool.h"
#include "CalXtalResponse/ICalTrigTool.h"
#include "../CalFailureMode/ICalFailureModeSvc.h"
#include "../CalCalib/IPrecalcCalibTool.h"

// GLAST
#include "geometry/Point.h"
#include "idents/VolumeIdentifier.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"
#include "Event/Digi/GltDigi.h"



// EXTLIB
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "TFile.h"
#include "TTree.h"

#include "CLHEP/Geometry/Transform3D.h"


// STD
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <set>
#include <cassert>
#include <map>
#include <memory>

using namespace CalUtil;
using namespace std;
using namespace CalibData;
using namespace idents;

//-- UTILITY METHODS --//

/// return relative diff (abs) between 2 floats avoid divide-by-zero errros
float rel_diff(float a, float b) {
  // safe to divide by a
  if (a!=0) return abs((a-b)/a);

  // fall back to divide by b
  if (b!=0) return abs((a-b)/b);

  // only possibility a==b==0
  // zero pct diff
  return 0;
}

/// \brief ensure that a threshold test is correct (given a certain margin for error)
/// \param test threshold level
/// \param input signal
/// \param measured test result
/// \param margin threhold to ignore if signal is close to margin. (in same units as threshold)

bool trig_test_margin(const float signal, const float thresh, const bool result, const double margin) {
  // return false if the trigger should have gone high, but it didn't
  if (signal >= thresh + margin && !result) return false;
  // return true if it should have gone low but fired anyway
  if (signal < thresh - margin && result) return false;

  // otherwise return true
  return true;
}

/// return true if C++ container contains given valueD
template <typename Iterator, typename Value>
bool contains(Iterator begin,
              Iterator end,
              const Value &val) {
  return std::find(begin, end, val) != end;
}
              
              

// Define the class here instead of in a header file: not needed anywhere but here!
//------------------------------------------------------------------------------
/** \class test_CalXtalRepsonse
    \brief Algorithm for unit_testing CalXtalResponse pkg functionality

    see execute() doc for test details
    
    \author Zach Fewtrell

*/

class test_CalXtalResponse : 
  public Algorithm {
public:
  test_CalXtalResponse(const string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

private: 
  void pos2Point(const XtalIdx xtalIdx, const float xtalPos, Point &pXtal);

  /// create single Mc hit & run it through full digi & recon
  /// check for expected output & accurate recon.
  StatusCode testSingleHit();

  /// for testing the outputs from a single xtal digi
  StatusCode verifyXtalDigi(const Event::CalDigi &calDigi);

  /// generate face signal & adc ped values from single xtal digi
  StatusCode preprocXtalDigi(const Event::CalDigi &calDigi);

  /// derive values from calib which are shared by several tests
  StatusCode preprocCalCalib();

  /// for verifying just FLE & FHE triggers from single xtal digi
  /// \param optional Glt object to check along w/ trigBits
  /// \param trigBits all trigger bits from one xtal to check
  StatusCode verifyXtalTrig(const Event::CalDigi &calDigi, 
                            CalArray<XtalDiode, bool> &trigBits,
                            Event::GltDigi const * const glt);

  /// \brief test calCalibSvc for current xtal & calib_src
  ///
  /// also load up calib values for use in analyzing higher level tests.
  StatusCode testCalCalibSvc();

  /// test all ULD, LAC, FLE/FHE threshold levels w/ series of hits
  StatusCode testTholds();

  /// test pedestal noise distrubution by reconning numerous identical deposits & checking the width.
  StatusCode testNoise();

  /// test HEX1 ADC channel saturation point
  StatusCode testSaturation();

  /// test sum of multiple
  StatusCode testMultiHit();


  /// fill mcIntegrating hit w/ given xtal position & energy
  void fillMcHit(const XtalIdx xtalIdx,
                 const float xtalPos,
                 const float meV,
                 Event::McIntegratingHit &hit);

  /// test CalFailureModeSvc
  StatusCode testCalFailureModeSvc();

  /// register McIntegratingHit col with event service.  remove & delete any existing hit col
  StatusCode registerMcHitCol(Event::McIntegratingHitCol *hitCol);

  /// zero out before each new test
  void newTest();

  /// store parameters for current test 
  /// so that i don't have to pass them from function to function
  class CurrentTest {
  public:
    CurrentTest(test_CalXtalResponse &test_cxr);

    /// zero out all values for new series of related tests
    void newTestSeq();

    /// currently tested xtal
    XtalIdx xtalIdx;
    /// current test energy
    float meV;
    /// current test deposit position in mm along crystal
    float xtalPos;
    /// current test number of hits.
    int nHits;

    /// string description of current test
    string testDesc;

    /// cal trig mode (BESTRNAGE / ALLRANGE) for current test
    CalXtalId::CalTrigMode trigMode;

    /// noise enabled for current test
    unsigned char noiseOn;
    /// zero suppression for current test
    unsigned char zeroSuppress;

    /// output lac threshold bits
    CalArray<FaceNum, bool> lacBits;
    /// output FLE/FHE trigger bits
    CalArray<XtalDiode, bool> trigBits;
    /// output ped subtracted adc values
    CalArray<XtalRng, float> adcPed;
    /// output estimated 'face signal' values
    CalArray<XtalRng, float> fsignl;

    /// output reconstructed energy
    float recEne;
    /// output reconstructed xtal position error
    float posDiff;

    /// current test calibration data source.
    enum CalibSrc {
      REAL_CALIB,
      IDEAL,
      N_CALIB_SRC
    };
    CalibSrc calibSrc;

    /// set calibration data source for current test
    void setCalibSrc(const CalibSrc src);

    /// get calibration data source for current test
    CalibSrc getCalibSrc() const {return calibSrc;}

    /// get current CalCalibSvc
    ICalCalibSvc  *getCalCalibSvc() const {return calCalibSvc;}
    /// get current XtalDigiTool
    IXtalDigiTool *getDigiTool() const {return digiTool;}
    /// get current CalSignalTool
    ICalSignalTool *getCalSignalTool() const {return (noiseOn) ? 
        signalToolNoise : signalTool;}
    /// get curren XtalRecTool
    IXtalRecTool  *getRecTool() const {return recTool;}
    /// get current CalTrigTool
    ICalTrigTool  *getTrigTool() const {return trigTool;}
    /// get current PrecalcCalibTool
    IPrecalcCalibTool *getPrecalcCalib() const {return precalcCalibTool;}

  private:    
    /// \brief used by constructor && reset f()'s
    ///
    /// clears out non-pointer primitives.
    /// \note does not set calibSrc
    void clear();

    /// cuyrrent CalCalibSvc
    ICalCalibSvc  *calCalibSvc;
    /// current XtalDigiTool (noise disabled)
    IXtalDigiTool *digiTool;
    /// current calsignal 
    ICalSignalTool *signalTool;
    ICalSignalTool *signalToolNoise;
    IXtalRecTool  *recTool;
    ICalTrigTool  *trigTool;
    IPrecalcCalibTool *precalcCalibTool;

    // use for reference to parent members.
    test_CalXtalResponse &parent;
  }; // CurrentTest

  /// point to current test info
  auto_ptr<CurrentTest> curTest;

  struct CurrentCalib {
    CurrentCalib() {clear();}

    // zero out all data members.
    void clear();

    CalArray<XtalRng, float>   ped;
    CalArray<XtalRng, float>   pedSig;
    CalArray<FaceNum, float>   lacThresh;
    CalArray<XtalRng, float>   uldThresh;

    CalArray<XtalDiode, float> trigMeV;
    CalArray<FaceNum, float> lacMeV;
    CalArray<XtalRng, float> uldMeV;
    /// mev value of 1 adc unit for each adc range.
    CalArray<XtalRng, float> mevPerADC;
  } curCalib;

  /// len of one Cal xtal
  float m_csiLength;

  /// the value of fLATObjects field, defining LAT towers 
  int m_eLATTowers; 
  
  /// the value of fTowerObject field, defining calorimeter module 
  int m_eTowerCAL;
  /// the value of fCellCmp field defining CsI crystal
  int m_eXtal;      
  /// number of geometric segments per Xtal
  int m_nCsISeg; 

  /// used for constants & conversion routines.
  IGlastDetSvc* m_detSvc;

  /// position estimation margin for single hit test (in mm)
  static const float m_singleHitPosMrgn;

  /// energy estimation margin for single hit test (in % of input ene)
  static const float m_singleHitEneMrgn;

  /// filename of output tuple.  No file created if set to default=""
  StringProperty m_tupleFilename;
  /// pointer to output tuple (TTree actually).  tuple is ignored
  /// if pointer is NULL
  TTree *m_tuple;
  /// pointer to tuple file.
  auto_ptr<TFile> m_tupleFile;

  /// by default, we only test a subset of xtals 
  /// by testing them all you can characterize a particular 
  /// calibration set.
  BooleanProperty m_testAllXtals;

  /// skip any tests w/ random generators, so we get same results every time.
  BooleanProperty m_skipNoise;

  /// 1000 per test would be perfect, but 100 per test is probably more practical.
  IntegerProperty m_nNoiseHits;

  set<TwrNum> m_testTwrs;
  set<LyrNum> m_testLyrs;
  set<ColNum> m_testCols;

  /// test ideal calibrations
  ICalCalibSvc *m_calCalibSvcIdeal;
  /// test ideal calibrations
  IXtalDigiTool *m_digiToolIdeal;
  /// test ideal calibrations
  IXtalDigiTool *m_digiToolIdealNoise;

  /// test ideal calibrations
  ICalSignalTool *m_signalToolIdeal;
  /// test ideal calibrations
  ICalSignalTool *m_signalToolIdealNoise;

  /// test ideal calibrations
  IXtalRecTool *m_recToolIdeal;
  /// test ideal calibrations
  ICalTrigTool *m_trigToolIdeal;
  /// test ideal calibrations
  IPrecalcCalibTool *m_precalcCalibToolIdeal;

  /// test  calibrations
  ICalCalibSvc *m_calCalibSvc;
  /// test  calibrations
  IXtalDigiTool *m_digiTool;
  /// test  calibrations
  ICalSignalTool *m_signalTool;
  /// test  calibrations
  ICalSignalTool *m_signalToolNoise;
  /// test  calibrations
  IXtalRecTool *m_recTool;
  /// test  calibrations
  ICalTrigTool *m_trigTool;
  /// test  calibrations
  IPrecalcCalibTool *m_precalcCalibTool;

  /// pointer to failure mode service
  ICalFailureModeSvc* m_calFailureModeSvc;
};

const float test_CalXtalResponse::m_singleHitEneMrgn = (float).01; // percent of input ene
const float test_CalXtalResponse::m_singleHitPosMrgn = 3; // (mm)

void test_CalXtalResponse::CurrentCalib::clear() {
  fill(ped.begin(), ped.end(), 0);
  fill(pedSig.begin(), pedSig.end(), 0);
  fill(lacThresh.begin(), lacThresh.end(), 0);
  fill(uldThresh.begin(), uldThresh.end(), 0);

  fill(trigMeV.begin(), trigMeV.end(), 0);
  fill(lacMeV.begin(), lacMeV.end(), 0);
  fill(uldMeV.begin(), uldMeV.end(), 0);
  fill(mevPerADC.begin(), mevPerADC.end(), 0);
}


test_CalXtalResponse::CurrentTest::CurrentTest(test_CalXtalResponse &cxr) :
  calCalibSvc(0),
  digiTool(0),
  signalTool(0),
  signalToolNoise(0),
  recTool(0),
  trigTool(0),
  precalcCalibTool(0),
  parent(cxr)
{
  clear();
}

void test_CalXtalResponse::CurrentTest::clear() {
  xtalIdx     = XtalIdx(0);
  meV         = 0;
  xtalPos     = 0;
  nHits       = 1;

  testDesc     = "";
  trigMode     = CalXtalId::BESTRANGE;
  noiseOn      = false;
  zeroSuppress = false;

  fill(lacBits.begin(), lacBits.end(), false);
  fill(trigBits.begin(), trigBits.end(), false);
  fill(adcPed.begin(), adcPed.end(), 0);
  fill(fsignl.begin(), fsignl.end(), 0);

  recEne = 0;
  posDiff = 0;
}

void test_CalXtalResponse::CurrentTest::newTestSeq() {
  clear();

  setCalibSrc(REAL_CALIB);
}

void test_CalXtalResponse::newTest() {
  curTest->getCalSignalTool()->newEvent();
  curTest->getCalSignalTool()->newEvent();
}

void test_CalXtalResponse::CurrentTest::setCalibSrc(const CalibSrc src) {
  calibSrc = src;
  switch (src) {
  case REAL_CALIB:
    calCalibSvc   = parent.m_calCalibSvc;
    digiTool      = parent.m_digiTool;
    signalTool      = parent.m_signalTool;
    signalToolNoise = parent.m_signalToolNoise;
    recTool       = parent.m_recTool;
    trigTool      = parent.m_trigTool;
    precalcCalibTool = parent.m_precalcCalibTool;
    break;
  case IDEAL:
    calCalibSvc   = parent.m_calCalibSvcIdeal;
    digiTool      = parent.m_digiToolIdeal;
    signalTool      = parent.m_signalToolIdeal;
    signalToolNoise = parent.m_signalToolIdealNoise;
    recTool       = parent.m_recToolIdeal;
    trigTool      = parent.m_trigToolIdeal;
    precalcCalibTool = parent.m_precalcCalibToolIdeal;
    break;
  default:
    assert(src >= 0 && src < N_CALIB_SRC);
  };
}

//------------------------------------------------------------------------

// necessary to define a Factory for this algorithm
// expect that the xxx_load.cxx file contains a call     
//     DLL_DECL_ALGORITHM( test_CalXtalResponse );

static const AlgFactory<test_CalXtalResponse>  Factory;
const IAlgFactory& test_CalXtalResponseFactory = Factory;

//------------------------------------------------------------------------
//! ctor
test_CalXtalResponse::test_CalXtalResponse(const string& name, ISvcLocator* pSvcLocator)
  : Algorithm(name, pSvcLocator),
    m_detSvc(0),
    m_tuple(0),
    m_calCalibSvcIdeal(0),
    m_digiToolIdeal(0),
    m_digiToolIdealNoise(0),
    m_signalToolIdeal(0),
    m_signalToolIdealNoise(0),
    m_recToolIdeal(0),
    m_trigToolIdeal(0),
    m_precalcCalibToolIdeal(0),
    m_calCalibSvc(0),
    m_digiTool(0),
    m_signalTool(0),
    m_signalToolNoise(0),
    m_recTool(0),
    m_trigTool(0),
    m_precalcCalibTool(0)

{
  declareProperty("tupleFilename",      m_tupleFilename      = "");
  declareProperty("testAllXtals",       m_testAllXtals       = false);
  declareProperty("nNoiseHits",         m_nNoiseHits         = 100);
  declareProperty("skipNoiseTests",     m_skipNoise          = false);

  curTest.reset(new CurrentTest(*this));
}

//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode test_CalXtalResponse::initialize(){
  StatusCode  sc;
  MsgStream msglog(msgSvc(), name()); 
  msglog << MSG::INFO << "initialize" << endreq;

  //-- RETRIEVE CONSTANTS --//
  // try to find the GlastDetSvc service
  sc = service("GlastDetSvc", m_detSvc);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't get GlastDetSvc " << endreq;
    return sc;
  }

  double tmp;
  // map containing pointers to integer constants to be read
  // with their symbolic names from xml file used as a key 
  typedef map<int*,string> PARAMAP;
  PARAMAP param;

  // INT CONSTANTS
  //     filling the map with information on constants to be read 
  param[&m_eTowerCAL]    = string("eTowerCAL");
  param[&m_eLATTowers]   = string("eLATTowers");
  param[&m_nCsISeg]      = string("nCsISeg");
  param[&m_eXtal]        = string("eXtal");

  // loop over all constants information contained in the map
  for(PARAMAP::iterator iter=param.begin(); iter!=param.end();iter++){
    //  retrieve constant
    if(!m_detSvc->getNumericConstByName((*iter).second, &tmp)) {
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::ERROR << " constant " <<(*iter).second
             << " not defined" << endreq;
      return StatusCode::FAILURE;
    } else *((*iter).first)= int(tmp); // store retrieved value 
  }

  // DOUBLE CONSTANTS
  sc = m_detSvc->getNumericConstByName("CsILength", &tmp);
  m_csiLength = tmp;
  if (sc.isFailure()) {
    msglog << MSG::ERROR << " constant CsILength not defined" << endreq;
    return sc;
  }


  //-- RETRIEVE TEST SERVICES --//
  sc = service("CalCalibSvcIdeal", m_calCalibSvcIdeal);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get CalCalibSvcIdeal." << endreq;
    return sc;
  }

  sc = service("CalCalibSvc", m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get CalCalibSvc." << endreq;
    return sc;
  }

  //-- RETRIEVE TEST TOOLS --//
  sc = toolSvc()->retrieveTool("XtalDigiTool", 
                               "XtalDigiToolIdeal",
                               m_digiToolIdeal,
                               this);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "XtalDigiToolIdeal" << endreq;
    return sc;
  }
  sc = toolSvc()->retrieveTool("CalSignalTool", 
                               "CalSignalToolIdeal",
                               m_signalToolIdeal,
                               0); // intended to be shared
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "CalSignalToolIdeal" << endreq;
    return sc;
  }
  sc = toolSvc()->retrieveTool("CalSignalTool", 
                               "CalSignalToolIdealNoise",
                               m_signalToolIdealNoise,
                               0); // intended to be shared
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "CalSignalToolIdeal" << endreq;
    return sc;
  }
  sc = toolSvc()->retrieveTool("XtalRecTool", 
                               "XtalRecToolIdeal",
                               m_recToolIdeal,
                               this);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "XtalRecToolIdeal" << endreq;
    return sc;
  }
  
  sc = toolSvc()->retrieveTool("CalTrigTool", 
                               "CalTrigToolIdeal",
                               m_trigToolIdeal,
                               this);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "CalTrigToolIdeal" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("PrecalcCalibTool", 
                               "PrecalcCalibToolIdeal",
                               m_precalcCalibToolIdeal,
                               0); // this tool is intended to be share by multiple users
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "PrecalcCalibToolIdeal" << endreq;
    return sc;
  }
 


  sc = toolSvc()->retrieveTool("CalTrigTool", 
                               "CalTrigTool",
                               m_trigTool,
                               this);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "CalTrigTool" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("PrecalcCalibTool", 
                               "PrecalcCalibTool",
                               m_precalcCalibTool,
                               0); // shared w/ other code
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "PrecalcCalibTool" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("XtalDigiTool", 
                               "XtalDigiTool",
                               m_digiTool,
                               this);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "XtalDigiTool" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("CalSignalTool", 
                               "CalSignalTool",
                               m_signalTool,
                               0); // intended to be shared
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "CalSignalTool" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("CalSignalTool", 
                               "CalSignalToolNoise",
                               m_signalToolNoise,
                               0); // intended to be shared
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "CalSignalTool" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("XtalRecTool", 
                               "XtalRecTool",
                               m_recTool,
                               this);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "XtalRecTool" << endreq;
    return sc;
  }


  // open optional tuple file
  if (m_tupleFilename.value().length() > 0 ) {
    m_tupleFile.reset(new TFile(m_tupleFilename.value().c_str(),"RECREATE","test_CalXtalResponse"));
    if (!m_tupleFile.get()) {
      msglog << MSG::ERROR << "Unable to create TTree object: " << m_tupleFilename << endreq;
      return StatusCode::FAILURE;
    }
    
    
    else {
      m_tuple = new TTree("test_CalXtalResponse","test_CalXtalResponse");
      if (!m_tuple) {
        msglog << MSG::ERROR << "Unable to create tuple" << endreq;
        return StatusCode::FAILURE;
      }

      //-- Add Branches to tree --//
      if (!m_tuple->Branch("xtalIdx", &curTest->xtalIdx, "xtalIdx/i") ||
          !m_tuple->Branch("meV", &curTest->meV, "meV/F") ||
          !m_tuple->Branch("xtalPos", &curTest->xtalPos, "xtalPos/F") ||
          !m_tuple->Branch("trigMode", &curTest->trigMode, "trigMode/b") ||
          !m_tuple->Branch("zeroSuppress", &curTest->zeroSuppress, "zeroSuppress/b") ||
          !m_tuple->Branch("noiseOn", &curTest->noiseOn, "noiseOn/b") ||
          !m_tuple->Branch("recEne", &curTest->recEne, "recEne/F") ||
          !m_tuple->Branch("posDiff", &curTest->posDiff, "posDiff/F") ||
          !m_tuple->Branch("calibSrc", &curTest->calibSrc, "calibSrc/b")
          ) {
        msglog << MSG::ERROR << "Couldn't create tuple branch" << endreq;
        return StatusCode::FAILURE;
      }
    }
  } // optional tuple

  //-- TEST XTAL SET (unless testAllXtals == true) --//
  // -- loop every test on each of these towers.
  m_testTwrs.insert(0);  // 1st
  m_testTwrs.insert(8);  // somewhere in middle, part of 8 tower partial lat.
  m_testTwrs.insert(15); // last

  // -- loop every test on each of these layers
  m_testLyrs.insert(LyrNum(0));  // 1st
  m_testLyrs.insert(LyrNum(4));  // somewhere in middle
  m_testLyrs.insert(LyrNum(7)); // las

  // -- loop every test on each of these columns
  m_testCols.insert(0);  // 1st
  m_testCols.insert(4);  // somewhere in middle
  m_testCols.insert(11);  // last

  sc = service("CalFailureModeSvc", m_calFailureModeSvc);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  Unable to find CalFailureMode service" << endreq;
    return sc;
  }
  
  return StatusCode::SUCCESS;
}

/** \brief Perform all CalXtalResponse tests in one enormous event

First run calFailureModeSvc

Now for the exhaustive test: 
- basically go through all of phase space, mc -> digi -> recon at each point & check recon results against mc deposit)
- noise is disabled in most tests to allow for precise recon / mc match.

Outer Loop - loop throuh set of select crystals spread out over different towers, layers and columns
 Inner Loop 1 - loop through calibration data sources (current ideal and real)
              - test CalCalibSvc for current crystal
  Inner Loop 2 - loop through full dynamic range of xtal (1,10,100,1000,10000 mev single xtal deposits
   Inner Loop 3 - loop through xtal longitudinal positions
                - test 4range readout
                - test bestrange
                - test zero suppression
                - enable noise (disabled in most tests)
   end loop 3
   - test noise widths by reconning numerous identical deposits.
   - test multiple hits along crystal.
   - test trigger / lac / uld thresholds
   - test for xtal saturation

 */
StatusCode test_CalXtalResponse::execute()
{
  MsgStream msglog(msgSvc(), name());
  StatusCode sc;
  
  //-- CalFailureModeSvc --//
  sc = testCalFailureModeSvc();
  if (sc.isFailure()) return sc;

  // -- test out following # of points evenly spaced along xtal
  const int nXtalPos(5); // -- gets absolute ends, 1 in ctr & 2 others
  
  // -- test out following energy levels (in a single deposit
  vector<float> testMeV;
  testMeV.push_back(0);
  testMeV.push_back(1); // below lac
  testMeV.push_back(1e1); // near muon peak
  testMeV.push_back(1e2); // near fle thresh
  testMeV.push_back(1e3);  // near fhe thresh, lex8 / hex1
  testMeV.push_back(1e4);

  for (XtalIdx xtalIdx; xtalIdx.isValid(); xtalIdx++) {
    //-- check if xtal is on test list --//
    if (!m_testAllXtals) {
      if (m_testTwrs.find(xtalIdx.getTwr()) == m_testTwrs.end()) continue;
      if (m_testLyrs.find(xtalIdx.getLyr()) == m_testLyrs.end()) continue;
      if (m_testCols.find(xtalIdx.getCol()) == m_testCols.end()) continue;
    }

    for (short calibSrc = 0; calibSrc < CurrentTest::N_CALIB_SRC; calibSrc++) {
      //////////////////////////
      //-- TEST CalCalibSvc --//
      //////////////////////////
      curTest->newTestSeq();
      curTest->xtalIdx = xtalIdx;
      curTest->setCalibSrc((CurrentTest::CalibSrc)calibSrc);
      sc = testCalCalibSvc();
      if (sc.isFailure()) return sc;
      


      //-- derive needed values from current xtal calib
      sc = preprocCalCalib();
      if (sc.isFailure()) return sc;

      for (vector<float>::const_iterator mevIt = testMeV.begin(); mevIt != testMeV.end(); mevIt++) {
        const float meV = *mevIt;
        for (int posIdx = 0; posIdx < nXtalPos; posIdx++) {

          //-- SHARED SETUP --//
          curTest->newTestSeq();
          curTest->xtalIdx = xtalIdx;
          curTest->meV     = meV;
          curTest->setCalibSrc((CurrentTest::CalibSrc)calibSrc);
          // longitudinal distance from ctr of xtal in mm
          curTest->xtalPos = m_csiLength*posIdx/(nXtalPos-1) - m_csiLength/2;


          ////////////////////////
          //-- 4RANGE READOUT --//
          ////////////////////////
          curTest->trigMode = CalXtalId::ALLRANGE;
          sc = testSingleHit();
          if (sc.isFailure()) return sc;
                  
          ///////////////////
          //-- BESTRANGE --//
          ///////////////////
          curTest->trigMode = CalXtalId::BESTRANGE;
          sc = testSingleHit();
          if (sc.isFailure()) return sc;     

          ///////////////////////
          //-- ZERO SUPPRESS --//
          ///////////////////////
          curTest->zeroSuppress = true;
          sc = testSingleHit();
          if (sc.isFailure()) return sc;
          curTest->zeroSuppress = false;
          
          if (!m_skipNoise) {
            //////////////////
            //-- NOISE ON --//
            //////////////////
            curTest->noiseOn = true;
            sc = testSingleHit();
            if (sc.isFailure()) return sc;
            curTest->noiseOn = false;
          }

        } // xtalPos loop

        if (!m_skipNoise) {
          //-- NOISE TEST --//
          sc = testNoise();
          if (sc.isFailure()) return sc;
        }

        //-- MULTI-HIT SUMMING TESTS --//
        sc = testMultiHit();
        if (sc.isFailure()) return sc;
          
      } // mev loop

      //-- TRIGGER THRESHOLD TESTS --//
      sc = testTholds();
      if (sc.isFailure()) return sc;

      sc = testSaturation();
      if (sc.isFailure()) return sc;

    } // calibSrc loop
  } // xtal loop  
  
  //-- TODO --//
  // EMPTY MC & DIGI TEST
  // DIODE DEPOSIT TEST
  // ALG TESTS
  // ZERO SUPPRESS TEST
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode test_CalXtalResponse::finalize(){
  MsgStream msglog(msgSvc(), name());

  // make sure optional tuple is closed out                                                        
  if (m_tupleFile.get()) {
    m_tupleFile->Write();
    m_tupleFile->Close(); // trees deleted                                                         
  }

  return StatusCode::SUCCESS;
}

/** \brief convert longitudinal xtalPos along cal xtal to 3d point in LAT geometry space

\param pXtal ouput position vector
\param xtalPos input longitudinal position in mm from center of xtal
*/
void test_CalXtalResponse::pos2Point(const XtalIdx xtalIdx, float xtalPos, Point &pXtal) {
  //-- CONSTRUCT GEOMETRY VECTORS FOR XTAL --//

  // create Volume Identifier for segments 0 & 11 of this crystal
  
  // volId info snagged from 
  // http://www.slac.stanford.edu/exp/glast/ground/software/geometry/docs/identifiers/geoId-RitzId.shtml
  idents::VolumeIdentifier segm0Id, segm11Id;

  const TwrNum twr = xtalIdx.getTwr();
  const LyrNum lyr = xtalIdx.getLyr();
  const ColNum col = xtalIdx.getCol();
  
  // init seg0 w/ info shared by both.
  segm0Id.append(m_eLATTowers);
  segm0Id.append(twr.getRow());
  segm0Id.append(twr.getCol());
  segm0Id.append(m_eTowerCAL);
  segm0Id.append(lyr.val());
  segm0Id.append(lyr.getDir().val()); 
  segm0Id.append(col.val());
  segm0Id.append(m_eXtal);

  // copy over shared info
  segm11Id = segm0Id;
  
  // init segment specific info.
  segm0Id.append(0); // segment Id
  segm11Id.append(m_nCsISeg-1); // set segment number for the last segment
  
  HepTransform3D transf;

  //get 3D transformation for segment 0 of this crystal
  m_detSvc->getTransform3DByID(segm0Id,&transf);
  //get position of the center of the segment 0
  Vector vect0 = transf.getTranslation();

  //get 3D transformation for the last segment of this crystal
  m_detSvc->getTransform3DByID(segm11Id,&transf);
  //get position of the center of the last segment
  Vector vect11 = transf.getTranslation();

  Point p0(0.,0.,0.);           
  // position of the crystal center
  Point pCenter = p0+(vect0+vect11)*0.5; 
  //normalized vector of the crystal direction 
  Vector dirXtal = (vect11-vect0)*m_nCsISeg/(m_nCsISeg-1);  

  // put 1D position info into 3D vector
  // 'xtalPos' is in units of xtal Length, convert to rel units (-1->1)
  xtalPos /= m_csiLength; 
  pXtal = pCenter+dirXtal*xtalPos;

}

/** \brief used in most tests.  full mc->digi->recon chain for single deposit
 */
StatusCode test_CalXtalResponse::testSingleHit() {
  StatusCode  sc;
  MsgStream msglog(msgSvc(), name()); 

  //-- init new test
  newTest();

  /////////////////////////////////////////////////
  // SECTION 0: GENERATE TEST DESCRIPTION STRING //
  /////////////////////////////////////////////////
  {
    ostringstream tmpStream;
    tmpStream << "TEST_SINGLE_HIT, "
              << curTest->xtalIdx.val() << ", "
              << curTest->meV           << ", "
              << curTest->xtalPos       << ", ";

    if (curTest->trigMode == CalXtalId::BESTRANGE)
      tmpStream << "BESTRANGE, ";
    else tmpStream << "ALLRANGE, ";

    if (curTest->zeroSuppress)
      tmpStream << "ZEROSUPPRESS, ";
    else tmpStream << "NOSUPPRESS, ";

    if (curTest->noiseOn) 
      tmpStream << "NOISEON, " ;
    else tmpStream << "NOISEOFF, ";

    switch (curTest->getCalibSrc()) {
    case CurrentTest::REAL_CALIB:
      tmpStream << "REAL_CALIB, ";
      break;
    case CurrentTest::IDEAL:
      tmpStream << "IDEAL_CALIB, ";
      break;
    default:
      short src = curTest->getCalibSrc();
      assert(src >= 0 && src < CurrentTest::N_CALIB_SRC);
    }

    curTest->testDesc = tmpStream.str();
  }

  //cout << "DEBUG: " << curTest->testDesc << endl;

  ////////////////////////////////////////////
  // SECTION 1: GENERATE MC INTEGRATING HIT //
  ////////////////////////////////////////////

  Event::McIntegratingHit *hit = new Event::McIntegratingHit();
  fillMcHit(curTest->xtalIdx, 
            curTest->xtalPos, 
            curTest->meV,
            *hit);

  //////////////////////////////
  //-- SECTION 2: XTAL DIGI --//
  //////////////////////////////
  // sum hit
  Event::McIntegratingHitCol *hitCol = new Event::McIntegratingHitCol();
  hitCol->push_back(hit);
  sc = registerMcHitCol(hitCol);
  if (sc.isFailure()) return sc;

  ICalSignalTool::XtalSignalMap xtalSignal;
  sc = curTest->getCalSignalTool()->getXtalSignalMap(curTest->xtalIdx, xtalSignal);
  if (sc.isFailure()) return sc;


  // digi outputs
  Event::CalDigi calDigi(curTest->trigMode, curTest->xtalIdx.getCalXtalId()); 
  sc = curTest->getDigiTool()->calculate(xtalSignal,
                                         calDigi,
                                         curTest->lacBits,
                                         curTest->zeroSuppress != 0);
  if (sc.isFailure()) return sc;
  
  
  //-- VERIFY XTAL DIGI OUTPUT --//
  sc = verifyXtalDigi(calDigi);
  if (sc.isFailure()) return sc;

  
  // recon outputs
  Event::CalXtalRecData recData(CalXtalId::BESTRANGE, curTest->xtalIdx.getCalXtalId());
  CalArray<FaceNum, bool> belowThresh;
  bool xtalBelowThresh;
  CalArray<FaceNum, bool> saturated;

  ///////////////////////////////////
  //-- SECTION 3: RUN XTAL RECON --//
  ///////////////////////////////////
  sc = curTest->getRecTool()->calculate(calDigi,
                                        recData,
                                        belowThresh,
                                        xtalBelowThresh,
                                        saturated,
                                        0);
  if (sc.isFailure()) return sc;

  // only process if ptr is non-null
  if (recData.getNReadouts() && !xtalBelowThresh) {
    Point testPoint;          // 3D point of test xtalPos
    pos2Point(curTest->xtalIdx, curTest->xtalPos, testPoint);
    curTest->posDiff = Vector(testPoint - recData.getPosition()).magnitude();
    curTest->recEne  = recData.getEnergy();
    const float eneDiff = rel_diff(curTest->meV, curTest->recEne);

    //--------- REPORT RESULT --------------------------//
    ostringstream tmp;
    tmp << curTest->testDesc 
        << setfill(' ') << setw(7) << setprecision(2) << curTest->posDiff << ", "
        << setfill(' ') << setw(7) << fixed << setprecision(2) << eneDiff << "%";
    curTest->testDesc = tmp.str();
    msglog << MSG::INFO << curTest->testDesc << endreq;

    //-- skip this particular test if noise is on. b/c answers will vary too much
    if (!curTest->noiseOn) {
      //--------- TEST MARGINS ---------------------------//
      if (curTest->posDiff > m_singleHitPosMrgn) {
        msglog << MSG::WARNING << "BAD POS, "
               << curTest->posDiff << ", "
               << curTest->xtalPos << ", "
               << curTest->testDesc << endreq;
        //return StatusCode::FAILURE;
      }
      
      if (eneDiff > m_singleHitEneMrgn) {
        msglog << MSG::WARNING << "BAD ENE, "
               << eneDiff << ", "
               << curTest->meV << ", "
               << curTest->testDesc << endreq;
        //return StatusCode::FAILURE;
      }
    }
  }

  if (m_tuple) 
    m_tuple->Fill();

  return StatusCode::SUCCESS;
}

StatusCode test_CalXtalResponse::verifyXtalDigi(const Event::CalDigi &calDigi) {
  StatusCode sc;
  MsgStream msglog(msgSvc(), name()); 
      
  //-- DIGI VERIFICTION: N READOUTS --//
  const short nRO = calDigi.getReadoutCol().size();

  //-- QUICK CHECK FOR ZERO SUPPRESSION (Exit early if need be) --//
  if (curTest->zeroSuppress) {
    // test that readout is only present if lacBits are true
    const bool anyLacTrue = contains(curTest->lacBits.begin(), 
                                     curTest->lacBits.end(), 
                                     true);
    if (anyLacTrue != (nRO != 0)) {
      msglog << MSG::ERROR << "BAD ZERO SUPPRESS, "
             << curTest->testDesc << endreq;
      return StatusCode::FAILURE;
    }

    // quit early if we (rightly so) have no readouts
    // as there is nothing left to test
    if (nRO == 0) 
      return StatusCode::SUCCESS;
  }

  const CalXtalId::CalTrigMode trigMode = calDigi.getMode();
  if ((trigMode == CalXtalId::BESTRANGE && nRO != 1) ||
      (trigMode == CalXtalId::ALLRANGE && nRO != 4)) {
    msglog << MSG::ERROR << "BAD N READOUTS, "
           << nRO << ", , "
           << curTest->testDesc << endreq;
    return StatusCode::FAILURE;
  }

  //-- PRE-PROCESS XTAL DIGI --//
  // generate ped subtracted adc & face Signal
  sc = preprocXtalDigi(calDigi);
  if (sc.isFailure()) return sc;

  // currently allways using 1st readout
  Event::CalDigi::CalXtalReadoutCol::const_iterator ro = 
    calDigi.getReadoutCol().begin();


  //-- DIGI VERIFICATION: --//
  for (FaceNum face; face.isValid(); face++) {
    const RngNum rng(ro->getRange(face));
    const XtalRng xRng(face,rng);

    //////////////
    // LAC test //
    //////////////
    if (!trig_test_margin(curTest->fsignl[xRng], 
                          curCalib.lacMeV[face], 
                          curTest->lacBits[face],
                          curCalib.mevPerADC[xRng]) 
        ) {
      msglog << MSG::ERROR << "BAD LAC BIT, "
             << curTest->fsignl[xRng] << ", "
             << curCalib.lacMeV[face] << ", "
             << curTest->testDesc << endreq;
      return StatusCode::FAILURE;
    }
    

    //////////////
    // ULD test //
    //////////////
    //-- ULD test1: val should be < uld limit for current range
    // ULD test1: val should be < uld limit for selected range (w/in 1 ADC UNIT)
    if (!trig_test_margin(curTest->adcPed[xRng], curCalib.uldThresh[xRng], false, 1)) {
      msglog << MSG::ERROR << "BAD ULD, "
             << curTest->adcPed[xRng] << ", "
             << curCalib.uldThresh[xRng] << ", "
             << curTest->testDesc << endreq;
      return StatusCode::FAILURE;
    }
 
    // ULD test2: val should be > uld limit for next range down
    if (rng > LEX8) {
      if (!curTest->noiseOn) //noise can mess up this test b/c it differs per channel
        if (!trig_test_margin(curTest->fsignl[xRng], 
                              curCalib.uldMeV[XtalRng(face, RngNum(rng.val()-1))], 
                              true, 
                              // xtra 1% margin for possible variation in faceSignal per channel
                              curCalib.mevPerADC[xRng] + .01*curTest->fsignl[xRng])) {
          msglog << MSG::ERROR << "BAD ULD, "
                 << curTest->fsignl[xRng] << ", "
                 << curCalib.uldMeV[XtalRng(face, RngNum(rng.val()-1))] << ", "
                 << curTest->testDesc << endreq;
          return StatusCode::FAILURE;
        }
    }
  } // face loop

  //-- CalTrigTool test --//
  CalArray<XtalDiode, bool> trigBitsTest;
  //-- clatrigtool pass 1: calDigi object, GLT off
  sc = curTest->getTrigTool()->calcXtalTrig(calDigi,
                                            trigBitsTest,
                                            NULL);
  if (sc.isFailure()) return sc;
  sc = verifyXtalTrig(calDigi, trigBitsTest, NULL);
  if (sc.isFailure()) return sc;

  //-- clatrigtool pass 2: GLT on
  Event::GltDigi glt;
  sc = curTest->getTrigTool()->calcXtalTrig(calDigi,
                                            trigBitsTest,
                                            &glt);
  if (sc.isFailure()) return sc;
  sc = verifyXtalTrig(calDigi, trigBitsTest, &glt);
  if (sc.isFailure()) return sc;

  //-- caltrigtool pass 3: this time w/ readout object & not calDigi object, GLT off --//
  fill(trigBitsTest.begin(), trigBitsTest.end(), false);
  sc = curTest->getTrigTool()->calcXtalTrig(curTest->xtalIdx,
                                            *ro,
                                            trigBitsTest,
                                            NULL);
  if (sc.isFailure()) return sc;
  sc = verifyXtalTrig(calDigi, trigBitsTest, NULL);
  if (sc.isFailure()) return sc;
                                           

  return StatusCode::SUCCESS;
}

/// check that trigger response matches energy deposits for single xtal.
StatusCode test_CalXtalResponse::verifyXtalTrig(const Event::CalDigi &calDigi,
                                                CalArray<XtalDiode, bool> &trigBits,
                                                Event::GltDigi const * const glt) {
  MsgStream msglog(msgSvc(), name()); 

  // currently allways using 1st readout
  Event::CalDigi::CalXtalReadoutCol::const_iterator ro = 
    calDigi.getReadoutCol().begin();


  for (FaceNum face; face.isValid(); face++) {
    const FaceIdx faceIdx(curTest->xtalIdx, face);

    const RngNum rng(ro->getRange(face));
    const XtalRng xRng(face,rng);
          
                                    
    //////////////
    // FLE test.//
    //////////////

    //-- any test based on estimated faceSignal can suffer due to noise since
    //   noise varies from channel to channel
    if (!curTest->noiseOn)
      if (!trig_test_margin(curTest->fsignl[xRng], 
                            curCalib.trigMeV[XtalDiode(face,LRG_DIODE)], 
                            trigBits[XtalDiode(face,LRG_DIODE)], 
                            // xtra 1% margin for possible variation in faceSignal per channel
                            curCalib.mevPerADC[xRng]+.01*curTest->fsignl[xRng])) {
        msglog << MSG::ERROR << "BAD FLE BIT, "
               << curTest->fsignl[xRng] << ", "
               << curCalib.trigMeV[XtalDiode(face,LRG_DIODE)] << ", "
               << curTest->testDesc << endreq;
        return StatusCode::FAILURE;
      }
        
    //////////////
    // FHE test.//
    //////////////
    //-- any test based on estimated faceSignal can suffer due to noise since
    //   noise varies from channel to channel
    if (!curTest->noiseOn)
      if (!trig_test_margin(curTest->fsignl[xRng], 
                            curCalib.trigMeV[XtalDiode(face,SM_DIODE)], 
                            trigBits[XtalDiode(face,SM_DIODE)], 
                            // xtra 1% margin for possible variation in faceSignal per channel
                            curCalib.mevPerADC[xRng]+.01*curTest->fsignl[xRng])) {
        msglog << MSG::WARNING            << "BAD FHE BIT, "
               << curTest->fsignl[xRng] << ", "
               << curCalib.trigMeV[XtalDiode(face,SM_DIODE)] << ", "
               << curTest->testDesc << endreq;
        //return StatusCode::FAILURE;
      }

    if (glt) {
      // test xtal FLE trigger
      if (trigBits[XtalDiode(face, LRG_DIODE)] !=
          glt->getCALLOtrigger(faceIdx.getCalXtalId())) {
        msglog << MSG::ERROR << "BAD GLT BIT, "
               << trigBits[XtalDiode(face, LRG_DIODE)] << ", "
               << glt->getCALLOtrigger(faceIdx.getCalXtalId()) << ", "
               << curTest->testDesc << endreq;
        return StatusCode::FAILURE;
      }
      
      // test xtal FHE trigger
      if (trigBits[XtalDiode(face, SM_DIODE)] !=
          glt->getCALHItrigger(faceIdx.getCalXtalId())) {
        msglog << MSG::ERROR << "BAD GLT BIT, "
               << trigBits[XtalDiode(face, SM_DIODE)] << ", "
               << glt->getCALHItrigger(faceIdx.getCalXtalId()) << ", "
               << curTest->testDesc << endreq;
        return StatusCode::FAILURE;
      }
    }
  }

  if (glt) {
    // test global FLE trigger
    if ((trigBits[XtalDiode(POS_FACE, LRG_DIODE)] ||
         trigBits[XtalDiode(NEG_FACE, LRG_DIODE)]) !=
        glt->getCALLOtrigger()) {
      msglog << MSG::ERROR << "BAD GLT BIT, "
             << (trigBits[XtalDiode(POS_FACE, LRG_DIODE)] ||
                 trigBits[XtalDiode(NEG_FACE, LRG_DIODE)]) << ", "
             << glt->getCALLOtrigger() << ", "
             << curTest->testDesc 
             << glt->getCALLOtrigger() << ", " << endreq;
      return StatusCode::FAILURE;
    }

    // test global FHE trigger
    if((trigBits[XtalDiode(POS_FACE, SM_DIODE)] ||
        trigBits[XtalDiode(NEG_FACE, SM_DIODE)]) != 
       glt->getCALHItrigger()) {
      msglog << MSG::ERROR << "BAD GLT BIT, "
             << (trigBits[XtalDiode(POS_FACE, SM_DIODE)] ||
                 trigBits[XtalDiode(NEG_FACE, SM_DIODE)]) << ", "
             << glt->getCALHItrigger() << ", "
             << curTest->testDesc << endreq;
      return StatusCode::FAILURE;
    }
  }

  return StatusCode::SUCCESS;
}


StatusCode test_CalXtalResponse::preprocXtalDigi(const Event::CalDigi &calDigi) {
  MsgStream msglog(msgSvc(), name());
  StatusCode sc;

  const int nReadouts = (curTest->trigMode == CalXtalId::BESTRANGE) ? 1 : 4;


  for (int nRO = 0; nRO < nReadouts; nRO++) {
    for (FaceNum face; face.isValid(); face++) {
      const float adc = calDigi.getAdc(nRO, face);
      const RngNum rng(calDigi.getRange(nRO, face));
      const RngIdx rngIdx(curTest->xtalIdx, face, rng);
      const XtalRng xRng(face,rng);
      
      curTest->adcPed[xRng] = adc - curCalib.ped[xRng];


      // calculate face signal from adc values
      sc = curTest->getCalCalibSvc()->evalFaceSignal(rngIdx, curTest->adcPed[xRng],
                                                     curTest->fsignl[xRng]);
      if (sc.isFailure()) return sc;

          
    }
  }

  return StatusCode::SUCCESS;
}


/// ensure that calibrations are all present and consistent for current crystal.
StatusCode test_CalXtalResponse::testCalCalibSvc() {
  StatusCode  sc;
  MsgStream msglog(msgSvc(), name()); 

  // get ready for new xtal's calibration
  curCalib.clear();

  /////////////////////////////////////////////////
  // SECTION 0: GENERATE TEST DESCRIPTION STRING //
  /////////////////////////////////////////////////
  {
    ostringstream tmpStream;
    tmpStream << "TEST_CALCALIBSVC, ";
    tmpStream << curTest->xtalIdx.val() << ", ";

    switch (curTest->getCalibSrc()) {
    case CurrentTest::REAL_CALIB:
      tmpStream << "REAL_CALIB, ";
      break;
    case CurrentTest::IDEAL:
      tmpStream << "IDEAL_CALIB, ";
      break;
    default:
      short src = curTest->getCalibSrc();
      assert(src >= 0 && src < CurrentTest::N_CALIB_SRC);
    }

    curTest->testDesc = tmpStream.str();
  }


  //-- MPD --//
  const CalMevPerDac *mpd;
  mpd = curTest->getCalCalibSvc()->getMPD(curTest->xtalIdx);
  if (!mpd) return StatusCode::FAILURE;

  //-- ASYM --//
  const CalAsym *asym;
  asym = curTest->getCalCalibSvc()->getAsym(curTest->xtalIdx);
  if (!asym) return StatusCode::FAILURE;

  CalArray<AsymType, const vector<ValSig> *> asymVals;

  asymVals[ASYM_LL] = asym->getBig();
  asymVals[ASYM_LS] = asym->getNSmallPBig();
  asymVals[ASYM_SL] = asym->getPSmallNBig();
  asymVals[ASYM_SS] = asym->getSmall();

  Xpos const*const xpos = curTest->getCalCalibSvc()->getAsymXpos();
  if (!xpos) return StatusCode::FAILURE;

  const vector<float> *xvals = xpos->getVals();
  
  if (xvals->size()   <= 0) {
    msglog << MSG::ERROR << "BAD ASYM, "
           << curTest->testDesc << endreq;
    return StatusCode::FAILURE;
  }

  for (unsigned i = 0; i < xvals->size(); i++) {
    for (AsymType asymType; asymType.isValid(); asymType++) {
      float testAsym, testPos;

      if (asymVals[asymType]->size() != xvals->size()) {
        msglog << MSG::ERROR << "BAD ASYM, "
               << curTest->testDesc << endreq;
        return StatusCode::FAILURE;
      }

      sc = curTest->getCalCalibSvc()->evalPos(curTest->xtalIdx, asymType, (*asymVals[asymType])[i].getVal(), testPos);
      if (sc.isFailure()) return sc;
      sc = curTest->getCalCalibSvc()->evalAsym(curTest->xtalIdx, asymType, (*xvals)[i], testAsym);
      if (sc.isFailure()) return sc;

      if (rel_diff((*asymVals[asymType])[i].getVal(), testAsym) > .01 ||
          rel_diff((*xvals)[i], testPos) > .01) {
        msglog << MSG::ERROR << "BAD ASYM, "
               << curTest->testDesc << endreq;
        return StatusCode::FAILURE;
      }
    
    } // asymtype
  }  // xVals

  //-- INTNONLIN && THOLDCI--//
  for (FaceNum face; face.isValid(); face++) {
    const FaceIdx faceIdx(curTest->xtalIdx, face);

    CalTholdCI const*const tholdCI = curTest->getCalCalibSvc()->getTholdCI(faceIdx);
    if (!tholdCI) return StatusCode::FAILURE;

    // load up new current calibrations
    curCalib.lacThresh[face]  = tholdCI->getLAC()->getVal();
        
    
    for (RngNum rng; rng.isValid(); rng++) {
      const XtalRng xRng(face,rng);
      const RngIdx rngIdx(curTest->xtalIdx, xRng);

      // load up new calibs
      curCalib.uldThresh[xRng] = tholdCI->getULD(rng.val())->getVal();
      
      const vector<float> *adcs = curTest->getCalCalibSvc()->getInlAdc(rngIdx);
      if (!adcs) return StatusCode::FAILURE;
      const vector<float> *cidacs = curTest->getCalCalibSvc()->getInlCIDAC(rngIdx);
      if (!cidacs) return StatusCode::FAILURE;

      if (cidacs->size() <= 0 || 
          adcs->size() > cidacs->size()) {
        msglog << MSG::ERROR << "BAD INL, "
               << curTest->testDesc << endreq;
        return StatusCode::FAILURE;
      }

      for (unsigned i = 0; i < adcs->size(); i++) {
        float testCIDAC;
        sc = curTest->getCalCalibSvc()->evalCIDAC(rngIdx, (*adcs)[i], testCIDAC);
        if (sc.isFailure()) return sc;

        if (rel_diff(testCIDAC, (*cidacs)[i]) > .01 && abs(testCIDAC - (*cidacs)[i] > .01)) {
          msglog << MSG::ERROR << "BAD INL, "
                 << curTest->testDesc << endreq;
          return StatusCode::FAILURE;
        }

        float testADC;
        sc = curTest->getCalCalibSvc()->evalADC(rngIdx, (*cidacs)[i], testADC);
        if (sc.isFailure()) return sc;

        if (rel_diff(testADC, (*adcs)[i]) > .01 && abs(testADC - (*adcs)[i] > .01)) {
          msglog << MSG::ERROR << "BAD INL, "
                 << curTest->testDesc << endreq;
          return StatusCode::FAILURE;
        }

      } // inl points

      Ped const*const ped = curTest->getCalCalibSvc()->getPed(rngIdx);
      if (!ped) return StatusCode::FAILURE;

      // load up new calibs
      curCalib.ped[xRng]    = ped->getAvr();
      curCalib.pedSig[xRng] = ped->getSig();
    }
  } // xRng
    
  return StatusCode::SUCCESS;
}


StatusCode test_CalXtalResponse::preprocCalCalib() {
  StatusCode  sc;
  MsgStream msglog(msgSvc(), name()); 
  
  for (FaceNum face; face.isValid(); face++) {
    // lacMeV
    const FaceIdx faceIdx(curTest->xtalIdx, face);
    sc = curTest->getCalCalibSvc()->evalFaceSignal(RngIdx(faceIdx,LEX8),
                                                   curCalib.lacThresh[face],
                                                   curCalib.lacMeV[face]);
    if (sc.isFailure()) return sc;

    // FLE & FHE MeV
    for (DiodeNum diode; diode.isValid(); diode++) {
      const XtalDiode xDiode(face, diode);
      const DiodeIdx diodeIdx(faceIdx, diode);
                                                   
      sc = curTest->getPrecalcCalib()->getTrigMeV(diodeIdx,
                                                  curCalib.trigMeV[xDiode]);
      if (sc.isFailure()) return sc;
    }
                                                  
    for (RngNum rng; rng.isValid(); rng++) {
      const RngIdx rngIdx(faceIdx,rng);
      const XtalRng xRng(face, rng);
      
      sc = curTest->getCalCalibSvc()->evalFaceSignal(rngIdx,
                                                     curCalib.uldThresh[xRng],
                                                     curCalib.uldMeV[xRng]);
      if (sc.isFailure()) return sc;
      
      // adcpermev - use last point in adc range
      curCalib.mevPerADC[xRng] = curCalib.uldMeV[xRng] / 
        curCalib.uldThresh[xRng];

    }
  }
  return StatusCode::SUCCESS;
}


/// simply run standard single hit tests at levels that will test each thold for on & off
/// testSingleHit() will check that all digital output in consistant including all triggers
/// & thresholds. (it does this by calling verifyXtalDigi()
StatusCode test_CalXtalResponse::testTholds() {
  StatusCode  sc;
  MsgStream msglog(msgSvc(), name()); 

  //-- all tests based on faceSignal, so use the center of
  //-- the xtal
  curTest->xtalPos = 0;

  //-- HEX1 SATURATION --//
  // run meV up to 95% of saturation for lowest saturating face
  // should work fine.
  curTest->meV = min(curCalib.uldMeV[XtalRng(POS_FACE,HEX1)],
                     curCalib.uldMeV[XtalRng(NEG_FACE,HEX1)])*.95;
  sc = testSingleHit();
  if (sc.isFailure()) return sc;
  

  for (FaceNum face; face.isValid(); face++) {
    //-- LAC --//
    curTest->meV = curCalib.lacMeV[face]*.95;
    sc = testSingleHit();
    if (sc.isFailure()) return sc;

    curTest->meV = curCalib.lacMeV[face]*1.05;
    sc = testSingleHit();
    if (sc.isFailure()) return sc;

    for (DiodeNum diode; diode.isValid(); diode++) {
      const XtalDiode xDiode(face,diode);
      //-- FLE/FHE --//
      curTest->meV = curCalib.trigMeV[xDiode]*.95;
      sc = testSingleHit();
      if (sc.isFailure()) return sc;

      curTest->meV = curCalib.trigMeV[xDiode]*1.05;
      sc = testSingleHit();
      if (sc.isFailure()) return sc;
      
    }

    //- HEX1 ULD is not a ULD, it's a saturation point.
    for (RngNum rng; rng < HEX1; rng++) {
      const XtalRng xRng(face,rng);
      //-- ULD --//
      curTest->meV = curCalib.uldMeV[xRng]*.95;
      sc = testSingleHit();
      if (sc.isFailure()) return sc;

      curTest->meV = curCalib.uldMeV[xRng]*1.05;
      sc = testSingleHit();
      if (sc.isFailure()) return sc;
      
    }
  }

  return StatusCode::SUCCESS;
}

void test_CalXtalResponse::fillMcHit(const XtalIdx xtalIdx,
                                     const float xtalPos,
                                     const float meV,
                                     Event::McIntegratingHit &hit) {

  //-- Create Volume Id (snagged from XtalRecTool::pos2Point()
  const TwrNum twr = xtalIdx.getTwr();
  const LyrNum lyr = xtalIdx.getLyr();
  const ColNum col = xtalIdx.getCol();
  idents::VolumeIdentifier segmId;
  segmId.append(m_eLATTowers);
  segmId.append(twr.getRow());
  segmId.append(twr.getCol());
  segmId.append(m_eTowerCAL);
  segmId.append(lyr.val());
  segmId.append(lyr.getDir().val()); 
  segmId.append(col.val());
  segmId.append(m_eXtal);

  //-- calculate xtal segment #
  short nSeg = (short)floor(((xtalPos + m_csiLength/2)  // range 0->CsILen(mm)
                             / m_csiLength)             // range 0->1
                            * m_nCsISeg);               // range 0->nCsISeg
            
  //-- make sure we clip the edges to segments 0 && (nSeg - 1)
  nSeg = max<short>(nSeg,0);
  nSeg = min<short>(nSeg,m_nCsISeg-1);

  //  apply seg # to volId
  segmId.append(nSeg);

  //-- calc 1st moment (distance from segment ctr) --//
  
  // position of segment center
  const float segCtr = (((float)nSeg+.5) // range 0.5 -> 11.5
                        / m_nCsISeg)     // range 0 -> 1
    * m_csiLength;   // range 0 -> csiLen
            
  // distance of hit from segment center
  const float xDiff = xtalPos + m_csiLength/2 - segCtr;

  const HepPoint3D mom1(xDiff,0,0);

  //-- Create McIntegratingHit
  hit.setVolumeID(segmId);
  // add single deposit w/ given mev & vector from center 
  // of segment.
  hit.addEnergyItem(meV, NULL, mom1);

}


StatusCode test_CalXtalResponse::testNoise() {
  StatusCode  sc;
  MsgStream msglog(msgSvc(), name()); 

  // test all ranges simultaneously
  curTest->trigMode = CalXtalId::ALLRANGE;
  curTest->nHits    = m_nNoiseHits;
  curTest->xtalPos  = 0;
  
  Event::McIntegratingHit hit;
  fillMcHit(curTest->xtalIdx, 
            curTest->xtalPos, 
            curTest->meV,
            hit);

  // -- first get 'noiseless' answer for comparison--//
  curTest->noiseOn = false;
  sc = testSingleHit();
  if (sc.isFailure()) return sc;
  CalArray<XtalRng, float> adcNoNoise;
  copy(curTest->adcPed.begin(), curTest->adcPed.end(),
       adcNoNoise.begin());
  
  //-- now run noise test --//
  curTest->noiseOn = true;
  
  // for calculating stdev & mean
  CalArray<XtalRng, double> adcSum;
  CalArray<XtalRng, double> adcSumSq;

  fill(adcSum.begin(), adcSum.end(), 0);
  fill(adcSumSq.begin(), adcSumSq.end(), 0);
  
  for (int nHit = 0; nHit < curTest->nHits; nHit++) {
    sc = testSingleHit();
    if (sc.isFailure()) return sc;
    for (XtalRng xRng; xRng.isValid(); xRng++) {
      //-- for each range, sum up adc values
      adcSum[xRng] += curTest->adcPed[xRng];
      adcSumSq[xRng] += pow(curTest->adcPed[xRng],2);      
    }
  }

  /////////////////////////////////////////////////
  // SECTION 0: GENERATE TEST DESCRIPTION STRING //
  /////////////////////////////////////////////////
  {
    ostringstream tmpStream;
    tmpStream << "TEST_NOISE, "
              << curTest->xtalIdx.val() << ", "
              << curTest->meV           << ", ";

    switch (curTest->getCalibSrc()) {
    case CurrentTest::REAL_CALIB:
      tmpStream << "REAL_CALIB, ";
      break;
    case CurrentTest::IDEAL:
      tmpStream << "IDEAL_CALIB, ";
      break;
    default:
      short src = curTest->getCalibSrc();
      assert(src >= 0 && src < CurrentTest::N_CALIB_SRC);
    }

    curTest->testDesc = tmpStream.str();
  }


  // calc & check 
  for (XtalRng xRng; xRng.isValid(); xRng++) {
    // do not look at ranges for which the ene level is too high
    // as saturation my bias results
    if (adcNoNoise[xRng] > curCalib.uldThresh[xRng]*.90) continue;
    // nor too low, as clipping may bias results
    if (adcNoNoise[xRng] < curCalib.pedSig[xRng]*3) continue;
    
    const float adcMean = adcSum[xRng]/curTest->nHits;
    const double meanSq  = adcSumSq[xRng]/curTest->nHits;
    const float adcRMS  = sqrt(meanSq - pow(adcMean,2));
    float diff = abs(adcRMS-curCalib.pedSig[xRng]);

    if (diff > 4*curCalib.pedSig[xRng]/sqrt((float)curTest->nHits)) {
      msglog << MSG::WARNING << "BAD NOISE SIGMA, "
             << adcRMS                << ", " 
             << curCalib.pedSig[xRng] << ", "
             << curTest->testDesc 
             << xRng.getFace().val() << ", "
             << xRng.getRng().val()  << ", " 
             << endreq;
      //return StatusCode::FAILURE;
    }

    const double thresh = 4*curCalib.pedSig[xRng]/sqrt((float)curTest->nHits);
    diff =  abs(adcMean-adcNoNoise[xRng]);
    if ( diff > thresh) {
      msglog << MSG::WARNING << "BAD NOISE MEAN, "
             << adcMean << ", "
             << adcNoNoise[xRng] << ", "
             << diff << ", "
             << thresh << ", "
             << curTest->testDesc 
             << xRng.getFace().val() << ", "
             << xRng.getRng().val()  << ", " 
             << endreq;
      //return StatusCode::FAILURE;
    }
  }

  return StatusCode::SUCCESS;
}

/// sum multiple mc hits across crystal -> digitize -> recon and check that recon matches.
StatusCode test_CalXtalResponse::testMultiHit() {
  StatusCode  sc;
  MsgStream msglog(msgSvc(), name()); 

  newTest();

  curTest->nHits   = 10;
  curTest->noiseOn = false;
  curTest->trigMode = CalXtalId::BESTRANGE;

  ////////////////////////
  // TEST B: ALONG XTAL //
  ////////////////////////
  // test several hits spread out along xtal
  curTest->xtalPos = 0; // average position will allways be at ctr 


  // to hold actual mc hits
  Event::McIntegratingHitCol *hitCol = new Event::McIntegratingHitCol();

  // create individual hits
  for (int nHit = 0; nHit < curTest->nHits; nHit++) {
    const float xtalPos = m_csiLength*nHit/(curTest->nHits-1) - m_csiLength/2;
    Event::McIntegratingHit *hit = new Event::McIntegratingHit();
    fillMcHit(curTest->xtalIdx, 
              xtalPos, 
              curTest->meV/curTest->nHits,
              *hit);
    hitCol->push_back(hit);    
  }

  /////////////////////////////////////////////////
  // SECTION 0: GENERATE TEST DESCRIPTION STRING //
  /////////////////////////////////////////////////
  {
    ostringstream tmpStream;
    tmpStream << "TEST_MULTIHIT_MULTIPOS, "
              << curTest->xtalIdx.val() << ", "
              << curTest->meV           << ", ";

    switch (curTest->getCalibSrc()) {
    case CurrentTest::REAL_CALIB:
      tmpStream << "REAL_CALIB, ";
      break;
    case CurrentTest::IDEAL:
      tmpStream << "IDEAL_CALIB, ";
      break;
    default:
      short src = curTest->getCalibSrc();
      assert(src >= 0 && src < CurrentTest::N_CALIB_SRC);
    }

    curTest->testDesc = tmpStream.str();
  }

  //-- XTAL DIGI --//
  
  // list of ptrs is needed for function call
  // build mc hit col
  sc = registerMcHitCol(hitCol);
  if (sc.isFailure())
    return sc;

  // retrieve signal values for given xtal
  ICalSignalTool::XtalSignalMap xtalSignal;
  sc = curTest->getCalSignalTool()->getXtalSignalMap(curTest->xtalIdx, xtalSignal);
  if (sc.isFailure()) return sc;

  Event::CalDigi calDigi(curTest->trigMode, curTest->xtalIdx.getCalXtalId()); 
  sc = curTest->getDigiTool()->calculate(xtalSignal,
                                         calDigi,
                                         curTest->lacBits,
                                         curTest->zeroSuppress != 0);
  if (sc.isFailure()) return sc;

  //-- VERIFY XTAL DIGI OUTPUT --//
  sc = verifyXtalDigi(calDigi);
  if (sc.isFailure()) return sc;

  // recon outputs
  Event::CalXtalRecData recData(CalXtalId::BESTRANGE, curTest->xtalIdx.getCalXtalId());
  CalArray<FaceNum, bool> belowThresh;
  bool xtalBelowThresh;
  CalArray<FaceNum, bool> saturated;

  ///////////////////////////////////
  //-- SECTION 3: RUN XTAL RECON --//
  ///////////////////////////////////
  sc = curTest->getRecTool()->calculate(calDigi,
                                        recData,
                                        belowThresh,
                                        xtalBelowThresh,
                                        saturated,
                                        0);
  if (sc.isFailure()) return sc;

  // only process if ptr is non-null
  if (recData.getNReadouts() && !xtalBelowThresh) {
    Point testPoint;          // 3D point of test xtalPos
    pos2Point(curTest->xtalIdx, curTest->xtalPos, testPoint);
    curTest->posDiff = Vector(testPoint - recData.getPosition()).magnitude();
    curTest->recEne  = recData.getEnergy();
    const float eneDiff = rel_diff(curTest->meV, curTest->recEne);
        
    //--------- REPORT RESULT --------------------------//
    ostringstream tmp;
    tmp << curTest->testDesc 
        << setfill(' ') << setw(7) << setprecision(2) << curTest->posDiff << ", "
        << setfill(' ') << setw(7) << fixed << setprecision(2) << eneDiff << "%";
    curTest->testDesc = tmp.str();
    msglog << MSG::INFO << curTest->testDesc << endreq;
        
    //--------- TEST MARGINS ---------------------------//
    if (curTest->posDiff > m_singleHitPosMrgn) {
      msglog << MSG::WARNING << "BAD POS, "
             << curTest->posDiff << ", "
             << curTest->xtalPos << ", "
             << curTest->testDesc << endreq;
      //return StatusCode::FAILURE;
    }
      
    if (eneDiff > m_singleHitEneMrgn) {
      msglog << MSG::WARNING << "BAD ENE, "
             << eneDiff << ", "
             << curTest->meV << ", "
             << curTest->testDesc << endreq;
      //return StatusCode::FAILURE;
    }
  }


  return StatusCode::SUCCESS;
}

/// simply put in too much energy for the xtal & the verifyXtalDigi() function
/// will ensure that the HEX1 adc values do not run too high.  Recon estimates
/// will be wron of course, so we'll just skip that part o the test.
StatusCode test_CalXtalResponse::testSaturation() {
  StatusCode  sc;
  MsgStream msglog(msgSvc(), name()); 

  //-- init for each trial
  newTest();

  /// inject 5% more energy than highest HEX1 saturation point for both faces.
  curTest->meV = 1.05*max(curCalib.uldMeV[XtalRng(POS_FACE, HEX1)],
                          curCalib.uldMeV[XtalRng(NEG_FACE, HEX1)]);
  // face signals are measured at center of xtal, so that's what we have to use.
  curTest->xtalPos = 0; 
  curTest->nHits = 1;
  curTest->noiseOn = false;
  

  /////////////////////////////////////////////////
  // SECTION 0: GENERATE TEST DESCRIPTION STRING //
  /////////////////////////////////////////////////
  {
    ostringstream tmpStream;
    tmpStream << "TEST_HEX1_SATURATION, ";
    tmpStream << curTest->xtalIdx.val() << ", ";
    tmpStream << curTest->meV << ", ";

    switch (curTest->getCalibSrc()) {
    case CurrentTest::REAL_CALIB:
      tmpStream << "REAL_CALIB, ";
      break;
    case CurrentTest::IDEAL:
      tmpStream << "IDEAL_CALIB, ";
      break;
    default:
      short src = curTest->getCalibSrc();
      assert(src >= 0 && src < CurrentTest::N_CALIB_SRC);
    }
    curTest->testDesc = tmpStream.str();
  }

  Event::McIntegratingHit *hit = new Event::McIntegratingHit();
  fillMcHit(curTest->xtalIdx, 
            curTest->xtalPos, 
            curTest->meV,
            *hit);


  //////////////////////////////
  //-- SECTION 2: XTAL DIGI --//
  //////////////////////////////
  // sum hit
  Event::McIntegratingHitCol *hitCol = new Event::McIntegratingHitCol();
  hitCol->push_back(hit);
  sc = registerMcHitCol(hitCol);
  if (sc.isFailure()) return sc;

  ICalSignalTool::XtalSignalMap xtalSignal;
  sc = curTest->getCalSignalTool()->getXtalSignalMap(curTest->xtalIdx, xtalSignal);
  if (sc.isFailure()) return sc;

  // digi outputs
  Event::CalDigi calDigi(curTest->trigMode, curTest->xtalIdx.getCalXtalId()); 
  sc = curTest->getDigiTool()->calculate(xtalSignal,
                                         calDigi,
                                         curTest->lacBits,
                                         curTest->zeroSuppress != 0);
  if (sc.isFailure()) return sc;
  
  
  //-- VERIFY XTAL DIGI OUTPUT --//
  sc = verifyXtalDigi(calDigi);
  if (sc.isFailure()) return sc;

  return StatusCode::SUCCESS;
}

StatusCode test_CalXtalResponse::testCalFailureModeSvc() {
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "testCalFailureModeSvc" << endreq;
  idents::CalXtalId id1(10,2,2);  // tower 10, layer 3(y)
  idents::CalXtalId id2(11,1,2);  // tower 11, layer 1(y)
  idents::CalXtalId id3(3,5,3);   // tower 3, layer 5(y)
  idents::CalXtalId id4(4,6,3);   // tower 4, layer 6(x)
  idents::CalXtalId id5(4,1,3);   // tower 4, layer 1(y)

  if (m_calFailureModeSvc == 0) return StatusCode::FAILURE;
  if (m_calFailureModeSvc->matchChannel(id1,idents::CalXtalId::NEG)) {
    log << MSG::INFO << "removed channel (10,2,2) NEG" << endreq;
  } else {
    log << MSG::ERROR << "failed to remove channel (10,2,2) NEG" << endreq;
    return StatusCode::FAILURE;
  }
  if (m_calFailureModeSvc->matchChannel(id2,idents::CalXtalId::NEG)) {
    log << MSG::INFO << "removed channel (11,1,2) NEG" << endreq;
  } else {
    log << MSG::ERROR << "failed to remove channel (11,1,2) NEG" << endreq;
    return StatusCode::FAILURE;
  }
  if (!m_calFailureModeSvc->matchChannel(id3,idents::CalXtalId::NEG)) {
    log << MSG::INFO << "left channel (3,5,3) NEG" << endreq;
  } else {
    log << MSG::ERROR << "erroneously removed channel (3,5,3) NEG" << endreq;
    return StatusCode::FAILURE;
  }
  if (m_calFailureModeSvc->matchChannel(id4,idents::CalXtalId::NEG)) {
    log << MSG::INFO << "removed channel (4,6,3) NEG" << endreq;
  } else {
    log << MSG::ERROR << "erroneously left channel (4,6,3) NEG" << endreq;
    return StatusCode::FAILURE;
  }
  if (!m_calFailureModeSvc->matchChannel(id5,idents::CalXtalId::NEG)) {
    log << MSG::INFO << "left channel (4,1,3) NEG" << endreq;
  } else {
    log << MSG::ERROR << "erroneously removed channel (4,1,3) NEG" << endreq;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

/// used to feed McHits to CalSignalSvc via TDS
StatusCode test_CalXtalResponse::registerMcHitCol(Event::McIntegratingHitCol *hitCol) {
  // check if object already exists
  DataObject *oldHitCol;
  StatusCode sc = eventSvc()->findObject(EventModel::MC::McIntegratingHitCol, oldHitCol);

  // unregister and destroy object if it currently resides in TDS
  if (oldHitCol != 0) {
    sc = eventSvc()->unregisterObject(EventModel::MC::McIntegratingHitCol);
    if (sc.isFailure())
      return sc;

    delete oldHitCol;
  }

  // create the TDS location for the McParticle Collection
  sc = eventSvc()->registerObject(EventModel::MC::McIntegratingHitCol, hitCol);
  if (sc.isFailure()) {
    MsgStream msglog(msgSvc(), name()); 
    msglog << "Failed to register McIntegratingHit" << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}
