// $Header$

// Include files
// Gaudi system includes

// LOCAL
#include "CalXtalResponse/ICalCalibSvc.h"
#include "CalXtalResponse/IXtalRecTool.h"
#include "CalXtalResponse/IXtalDigiTool.h"
#include "CalXtalResponse/ICalTrigTool.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

// EXTLIB
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"
#include "TTree.h"

// STD
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <set>
#include <cassert>
#include <cstring>

using namespace CalUtil;
using namespace std;
using namespace CalibData;

/// return % diff (abs) between 2 floats avoid divide-by-zero errros
float pct_diff(float a, float b) {
  // safe to divide by a
  if (a!=0) return abs(a-b)/a;

  // fall back to divide by b
  if (b!=0) return abs(a-b)/b;

  // only possibility a==b==0
  // zero pct diff
  return 0;
}

/// \brief ensure that a threshold test is correct (given a certain margin for error)
/// \param test threshold level
/// \param input signal
/// \param measured test result
/// \param margin threhold to ignore if signal is close to margin. (in same units as threshold)

bool trig_test_margin(float signal, float thresh, bool result, double margin) {
  // return false if the trigger should have gone high, but it didn't
  if (signal >= thresh + margin && !result) return false;
  // return true if it should have gone low but fired anyway
  if (signal < thresh - margin && result) return false;

  // otherwise return true
  return true;
}

/// \brief ensure that a threshold test is correct (given a certain relative margin for error)
/// \param test threshold level
/// \param input signal
/// \param measured test result
/// \param margin threhold to ignore if signal is close to margin. (as a fraction of the threhold

bool trig_test_pct_margin(float signal, float thresh, bool result, double margin) {
  // return false if the trigger should have gone high, but it didn't
  if (signal >= thresh + margin * thresh && !result) return false;
  // return true if it should have gone low but fired anyway
  if (signal < thresh - thresh * margin && result) return false;

  // otherwise return true
  return true;
}

// Define the class here instead of in a header file: not needed anywhere but here!
//------------------------------------------------------------------------------
/** 
    A simple algorithm.

*/

/// test ideal calibrations
ICalCalibSvc *m_calCalibSvcIdeal=0;
/// test ideal calibrations
IXtalDigiTool *m_digiToolIdeal=0;
IXtalDigiTool *m_digiToolIdealNoise=0;

/// test ideal calibrations
IXtalRecTool *m_recToolIdeal=0;
/// test ideal calibrations
ICalTrigTool *m_trigToolIdeal=0;


/// test real calibrations
ICalCalibSvc *m_calCalibSvcReal=0;
/// test real calibrations
IXtalDigiTool *m_digiToolReal=0;
IXtalDigiTool *m_digiToolRealNoise=0;
/// test real calibrations
IXtalRecTool *m_recToolReal=0;
/// test real calibrations
ICalTrigTool *m_trigToolReal=0;

/// test real calibrations
ICalCalibSvc *m_calCalibSvc8Tower=0;
/// test 8Tower calibrations
IXtalDigiTool *m_digiTool8Tower=0;
IXtalDigiTool *m_digiTool8TowerNoise=0;
/// test 8Tower calibrations
IXtalRecTool *m_recTool8Tower=0;
/// test 8Tower calibrations
ICalTrigTool *m_trigTool8Tower=0;

class test_CalXtalResponse : public Algorithm {
public:
  test_CalXtalResponse(const string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

private: 
  /// count each invokation of execute() function.
  int m_count;

  /// length of one Cal xtal
  float m_CsILength;

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
  


  void pos2Point(XtalIdx xtalIdx, float xtalPos, Point &pXtal);
  StatusCode hex8_to_hex1(FaceIdx faceIdx, float h8adc, float &h1adc, 
                          ICalCalibSvc &calCalibSvc);
  StatusCode lex8_to_lex1(FaceIdx faceIdx, float l8adc, float &l1adc, 
                          ICalCalibSvc &calCaibSvc);


  /// create single Mc hit & run it through full digi & recon
  /// check for expected output & accurate recon.
  StatusCode testSingleHit();

  /// for testing the outputs from a single xtal digi
  StatusCode verifyXtalDigi(const Event::CalDigi &calDigi, Event::GltDigi *glt);

  /// generate face signal & adc ped values from single xtal digi
  StatusCode preprocXtalDigi(const Event::CalDigi &calDigi);

  /// derive values from calib which are shared by several tests
  StatusCode preprocCalCalib();

  /// for verifying just FLE & FHE triggers from single xtal digi
  /// \param optional Glt object to check along w/ trigBits
  /// \param trigBits all trigger bits from one xtal to check
  StatusCode verifyXtalTrig(const Event::CalDigi &calDigi, 
                            CalArray<XtalDiode, bool> &trigBits,
                            Event::GltDigi *glt);

  /// test calCalibSvc for current xtal & calib_src
  StatusCode testCalCalibSvc();

  /// test all ULD, LAC, FLE/FHE threshold levels w/ series of hits
  StatusCode testTholds();

  /// test pedestal noise distrubution
  StatusCode testNoise();

  /// test HEX1 ADC channel saturation point
  StatusCode testSaturation();

  /// test sum of multiple
  StatusCode testMultiHit();

  /// fill mcIntegrating hit w/ given xtal position & energy
  void fillMcHit(XtalIdx xtalIdx,
                 float xtalPos,
                 float meV,
                 Event::McIntegratingHit &hit);


  /// position estimation margin for single hit test (in mm)
  static const float m_singleHitPosMrgn;

  /// energy estimation margin for single hit test (in % of input ene)
  static const float m_singleHitEneMrgn;

  /// set of tower bay numbers included in 8 tower geometry
  set<TwrNum> m_8towerSet;

private:
  /// store parameters for current test 
  /// so that i don't have to pass them from function to function
  struct CurrentTest {
    void CurentTest() {clear();}

    /// zero out all values for new test
    void clear();

    XtalIdx xtalIdx;
    float meV;
    float xtalPos;
    int nHits;

    string testDesc;

    CalXtalId::CalTrigMode trigMode;

    unsigned char tupleOn;
    unsigned char gltOn;
    unsigned char noiseOn;

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

    enum CalibSrc {
      FLIGHT_LAT,
      PARTIAL_LAT,
      IDEAL,
      N_CALIB_SRC
    };

    void setCalibSrc(CalibSrc src);
    CalibSrc getCalibSrc() const {return calibSrc;}
    ICalCalibSvc  *getCalCalibSvc() const {return calCalibSvc;}
    IXtalDigiTool *getDigiTool() const {return (noiseOn) ? 
                                          digiToolNoise : digiTool;}
    IXtalRecTool  *getRecTool() const {return recTool;}
    ICalTrigTool  *getTrigTool() const {return trigTool;}

    CalibSrc calibSrc;
  private:    
    ICalCalibSvc  *calCalibSvc;
    IXtalDigiTool *digiTool;
    IXtalDigiTool *digiToolNoise;
    IXtalRecTool  *recTool;
    ICalTrigTool  *trigTool;

  } curTest;

  struct CurrentCalib {
    CurrentCalib() {clear();}
    void clear() {memset(this,0,sizeof(this));}
    CalArray<XtalRng, float>   ped;
    CalArray<XtalRng, float>   pedSig;
    CalArray<XtalDiode, float> trigThresh;
    CalArray<FaceNum, float>   lacThresh;
    CalArray<XtalRng, float>   uldThresh;

    CalArray<XtalDiode, float> trigMeV;
    CalArray<FaceNum, float> lacMeV;
    CalArray<XtalRng, float> uldMeV;
    /// mev value of 1 adc unit for each adc range.
    CalArray<XtalRng, float> mevPerADC;
    /// pedestal noise in MeV
    CalArray<XtalRng, float> pedNoise;
  } curCalib;

  /// filename of output tuple.  No file created if set to default=""
  StringProperty m_tupleFilename;
  /// pointer to output tuple (TTree actually).  tuple is ignored
  /// if pointer is NULL
  TTree *m_tuple;
  /// pointer to XtalDigiToolTuple file.
  TFile *m_tupleFile;

  /// by default, we only test a subset of xtals 
  /// by testing them all you can characterize a particular 
  /// calibration set.
  BooleanProperty m_testAllXtals;

  /// 1000 per test would be perfect, but 100 per test is probably more practical.
  IntegerProperty m_nNoiseHits;

  set<TwrNum> m_testTwrs;
  set<LyrNum> m_testLyrs;
  set<ColNum> m_testCols;
};

const float test_CalXtalResponse::m_singleHitEneMrgn = (float).01; // percent of input ene
const float test_CalXtalResponse::m_singleHitPosMrgn = 3; // (mm)

void test_CalXtalResponse::CurrentTest::clear() {
  xtalIdx     = XtalIdx(0);
  meV         = 0;
  xtalPos     = 0;
  nHits       = 1;

  testDesc    = "";
  trigMode    = CalXtalId::BESTRANGE;
  tupleOn     = true;
  gltOn       = true;
  noiseOn     = false;

  lacBits.fill(false);
  trigBits.fill(false);
  adcPed.fill(0);
  fsignl.fill(0);

  recEne = 0;
  posDiff = 0;

  setCalibSrc(FLIGHT_LAT);
}

void test_CalXtalResponse::CurrentTest::setCalibSrc(CalibSrc src) {
  calibSrc = src;
  switch (src) {
  case FLIGHT_LAT:
    calCalibSvc   = m_calCalibSvcReal;
    digiTool      = m_digiToolReal;
    digiToolNoise = m_digiToolRealNoise;
    recTool       = m_recToolReal;
    trigTool      = m_trigToolReal;
    break;
  case PARTIAL_LAT:
    calCalibSvc   = m_calCalibSvc8Tower;
    digiTool      = m_digiTool8Tower;
    digiToolNoise = m_digiTool8TowerNoise;
    recTool       = m_recTool8Tower;
    trigTool      = m_trigTool8Tower;
    break;
  case IDEAL:
    calCalibSvc   = m_calCalibSvcIdeal;
    digiTool      = m_digiToolIdeal;
    digiToolNoise = m_digiToolIdealNoise;
    recTool       = m_recToolIdeal;
    trigTool      = m_trigToolIdeal;
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
    m_count(0),
    m_detSvc(0),
    m_tuple(0),
    m_tupleFile(0)
{
  m_8towerSet.insert(0);
  m_8towerSet.insert(1);
  m_8towerSet.insert(4);
  m_8towerSet.insert(5);
  m_8towerSet.insert(8);
  m_8towerSet.insert(9);
  m_8towerSet.insert(12);
  m_8towerSet.insert(13);

  declareProperty("tupleFilename",      m_tupleFilename      = "");
  declareProperty("testAllXtals",       m_testAllXtals       = false);
  declareProperty("nNoiseHits",         m_nNoiseHits         = 100);

}

//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode test_CalXtalResponse::initialize(){
  StatusCode  sc;
  MsgStream msglog(msgSvc(), name()); 
  msglog << MSG::INFO << "initialize" << endreq;

  //-- RETRIEVE CONSTANTS --//
  // try to find the GlastDevSvc service
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
             <<" not defined" << endreq;
      return StatusCode::FAILURE;
    } else *((*iter).first)= int(tmp); // store retrieved value 
  }

  // DOUBLE CONSTANTS
  sc = m_detSvc->getNumericConstByName("CsILength", &tmp);
  m_CsILength = tmp;
  if (sc.isFailure()) {
    msglog << MSG::ERROR << " constant CsILength not defined" << endreq;
    return StatusCode::FAILURE;
  }


  //-- RETRIEVE TEST SERVICES --//
  sc = service("CalCalibSvcIdeal", m_calCalibSvcIdeal);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get CalCalibSvcIdeal." << endreq;
    return sc;
  }

  sc = service("CalCalibSvcReal", m_calCalibSvcReal);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get CalCalibSvcReal." << endreq;
    return sc;
  }
  
  sc = service("CalCalibSvc8Tower", m_calCalibSvc8Tower);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get CalCalibSvc8Tower." << endreq;
    return sc;
  }

  //-- RETRIEVE TEST TOOLS --//
  sc = toolSvc()->retrieveTool("XtalDigiTool", 
                               "XtalDigiToolIdeal",
                               m_digiToolIdeal);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "XtalDigiToolIdeal" << endreq;
    return sc;
  }
  sc = toolSvc()->retrieveTool("XtalDigiTool", 
                               "XtalDigiToolIdealNoise",
                               m_digiToolIdealNoise);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "XtalDigiToolIdeal" << endreq;
    return sc;
  }
  sc = toolSvc()->retrieveTool("XtalRecTool", 
                               "XtalRecToolIdeal",
                               m_recToolIdeal);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "XtalRecToolIdeal" << endreq;
    return sc;
  }
  
  sc = toolSvc()->retrieveTool("CalTrigTool", 
                               "CalTrigToolIdeal",
                               m_trigToolIdeal);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "CalTrigToolIdeal" << endreq;
    return sc;
  }
 


  sc = toolSvc()->retrieveTool("CalTrigTool", 
                               "CalTrigToolReal",
                               m_trigToolReal);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "CalTrigToolReal" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("XtalDigiTool", 
                               "XtalDigiToolReal",
                               m_digiToolReal);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "XtalDigiToolReal" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("XtalDigiTool", 
                               "XtalDigiToolRealNoise",
                               m_digiToolRealNoise);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "XtalDigiToolReal" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("XtalRecTool", 
                               "XtalRecToolReal",
                               m_recToolReal);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "XtalRecToolReal" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("CalTrigTool", 
                               "CalTrigTool8Tower",
                               m_trigTool8Tower);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "CalTrigTool8Tower" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("XtalDigiTool", 
                               "XtalDigiTool8Tower",
                               m_digiTool8Tower);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "XtalDigiTool8Tower" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("XtalDigiTool", 
                               "XtalDigiTool8TowerNoise",
                               m_digiTool8TowerNoise);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "XtalDigiTool8Tower" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("XtalRecTool", 
                               "XtalRecTool8Tower",
                               m_recTool8Tower);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "XtalRecTool8Tower" << endreq;
    return sc;
  }


  // open optional tuple file
  if (m_tupleFilename.value().length() > 0 ) {
    m_tupleFile = new TFile(m_tupleFilename.value().c_str(),"RECREATE","test_CalXtalResponse");
    if (!m_tupleFile)
      // allow to continue w/out tuple file as it is not a _real_ failure
      msglog << MSG::ERROR << "Unable to create TTree object: " << m_tupleFilename << endreq;
    
    else {
      m_tuple = new TTree("test_CalXtalResponse","test_CalXtalResponse");
      if (!m_tuple) {
        msglog << MSG::ERROR << "Unable to create tuple" << endreq;
        return StatusCode::FAILURE;
      }

      //-- Add Branches to tree --//
      if (!m_tuple->Branch("xtalIdx", &curTest.xtalIdx, "xtalIdx/i") ||
          !m_tuple->Branch("meV", &curTest.meV, "meV/F") ||
          !m_tuple->Branch("xtalPos", &curTest.xtalPos, "xtalPos/F") ||
          !m_tuple->Branch("trigMode", &curTest.trigMode, "trigMode/b") ||
          !m_tuple->Branch("tupleOn", &curTest.tupleOn, "tupleOn/b") ||
          !m_tuple->Branch("gltOn", &curTest.gltOn, "gltOn/b") ||
          !m_tuple->Branch("noiseOn", &curTest.noiseOn, "noiseOn/b") ||
          !m_tuple->Branch("recEne", &curTest.recEne, "recEne/F") ||
          !m_tuple->Branch("posDiff", &curTest.posDiff, "posDiff/F") ||
          !m_tuple->Branch("calibSrc", &curTest.calibSrc, "calibSrc/b")
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
  m_testLyrs.insert(0);  // 1st
  m_testLyrs.insert(4);  // somewhere in middle
  m_testLyrs.insert(7); // las

  // -- loop every test on each of these columns
  m_testCols.insert(0);  // 1st
  m_testCols.insert(4);  // somewhere in middle
  m_testCols.insert(11);  // last

  
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------
//! process an event
StatusCode test_CalXtalResponse::execute()
{
  MsgStream msglog(msgSvc(), name());
  msglog << MSG::INFO << "executing " << ++m_count << " time" << endreq;
  StatusCode sc;
  

  // -- test out following # of points evenly spaced along xtal
  int nXtalPos(5); // -- gets absolute ends, 1 in ctr & 2 others
  
  // -- test out following energy levels (in a single deposit
  vector<float> testMeV;
  testMeV.push_back(1);
  testMeV.push_back(1e1);
  testMeV.push_back(1e2);
  testMeV.push_back(1e3);  // 1 GeV // currently running in muon mode, so we'll max out at 1gev per xtal

  // > 70 GeV should saturate xtal.... will do this in later test.
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
      curTest.clear();
      curTest.xtalIdx = xtalIdx;
      curTest.setCalibSrc((CurrentTest::CalibSrc)calibSrc);
      sc = testCalCalibSvc();
      //-- PARTIAL LAT TEST --//
      if (curTest.getCalibSrc() == curTest.PARTIAL_LAT) {
        // SHOULD _NOT_ FAIL
        if (m_8towerSet.find(curTest.xtalIdx.getTwr()) != m_8towerSet.end()) {
          if (sc.isFailure()) return sc;
        } else {// _SHOULD_ FAIL
          if (!sc.isFailure()) {
            msglog << MSG::ERROR << "TESTFAIL, BAD PARTIAL LAT, "
                   << curTest.xtalIdx.getTwr() << ", , "
                   << curTest.testDesc << endreq;
            return StatusCode::FAILURE;
          }

          // quietly skip rest of testing on this xtal as it is not installed in LAT.
          continue;
        }
      }
      else //-- SHOULD BE NO FAILURES FOR 16 TOWER LAT --//
        if (sc.isFailure()) return sc;

      //-- derive needed values from current xtal calib
      sc = preprocCalCalib();
      if (sc.isFailure()) return sc;

      for (vector<float>::const_iterator mevIt = testMeV.begin(); mevIt != testMeV.end(); mevIt++) {
        float meV = *mevIt;
        for (int posIdx = 0; posIdx < nXtalPos; posIdx++) {
          // longitudinal distance from ctr of xtal in mm
          curTest.xtalPos = m_CsILength*posIdx/(nXtalPos-1) - m_CsILength/2;


          //-- SHARED SETUP --//
          curTest.clear();
          curTest.xtalIdx = xtalIdx;
          curTest.meV     = meV;
          curTest.setCalibSrc((CurrentTest::CalibSrc)calibSrc);

          ////////////////////////
          //-- 4RANGE READOUT --//
          ////////////////////////
          curTest.trigMode= CalXtalId::ALLRANGE;
          sc = testSingleHit();
          if (sc.isFailure()) return sc;
                  
          ///////////////////
          //-- BESTRANGE --//
          ///////////////////
          curTest.trigMode= CalXtalId::BESTRANGE;
          sc = testSingleHit();
          if (sc.isFailure()) return sc;     

          //////////////////////
          //-- CalTuple OFF --//
          //////////////////////
          curTest.tupleOn = false;
          sc = testSingleHit();
          if (sc.isFailure()) return sc;
          curTest.tupleOn = true;

          /////////////////
          //-- GLT OFF --//
          /////////////////
          curTest.gltOn = false;
          sc = testSingleHit();
          if (sc.isFailure()) return sc;
          curTest.gltOn = true;
          
          //////////////////
          //-- NOISE ON --//
          //////////////////
          curTest.noiseOn = true;
          sc = testSingleHit();
          if (sc.isFailure()) return sc;
          curTest.noiseOn = false;

        } // xtalPos loop

        //-- NOISE TEST --//
        sc = testNoise();
        if (sc.isFailure()) return sc;

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
  // ZERO ENE TEST
  // DIODE DEPOSIT TEST
  // CALXTALRECALG TEST
  // CALDIGIALG TEST
  // CALXTALRESP W/ OUT LAC CUT
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode test_CalXtalResponse::finalize(){
  MsgStream msglog(msgSvc(), name());
  msglog << MSG::INFO << "finalize after " << m_count << " calls." << endreq;

  // make sure optional tuple is closed out                                                        
  if (m_tupleFile) {
    m_tupleFile->Write();
    m_tupleFile->Close(); // trees deleted                                                         
    delete m_tupleFile;
    m_tupleFile = 0;
  }


  return StatusCode::SUCCESS;
}

/** \brief convert longitudinal xtalPos along cal xtal to 3d point in LAT geometry space

\param pXtal ouput position vector
\param xtalPos input longitudinal position in mm from center of xtal
*/
void test_CalXtalResponse::pos2Point(XtalIdx xtalIdx, float xtalPos, Point &pXtal) {
  //-- CONSTRUCT GEOMETRY VECTORS FOR XTAL --//

  // create Volume Identifier for segments 0 & 11 of this crystal
  
  // volId info snagged from 
  // http://www.slac.stanford.edu/exp/glast/ground/software/geometry/docs/identifiers/geoId-RitzId.shtml
  idents::VolumeIdentifier segm0Id, segm11Id;

  TwrNum twr = xtalIdx.getTwr();
  LyrNum lyr = xtalIdx.getLyr();
  ColNum col = xtalIdx.getCol();
  
  // init seg0 w/ info shared by both.
  segm0Id.append(m_eLATTowers);
  segm0Id.append(twr.getRow());
  segm0Id.append(twr.getCol());
  segm0Id.append(m_eTowerCAL);
  segm0Id.append(lyr);
  segm0Id.append(lyr.getDir()); 
  segm0Id.append(col);
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
  xtalPos /= m_CsILength; 
  pXtal = pCenter+dirXtal*xtalPos;

}

StatusCode test_CalXtalResponse::testSingleHit() {
  StatusCode  sc;
  MsgStream msglog(msgSvc(), name()); 

  TwrNum twr = curTest.xtalIdx.getTwr();
  LyrNum lyr = curTest.xtalIdx.getLyr();
  ColNum col = curTest.xtalIdx.getCol();

  /////////////////////////////////////////////////
  // SECTION 0: GENERATE TEST DESCRIPTION STRING //
  /////////////////////////////////////////////////
  {
    ostringstream tmpStream;
    tmpStream << "TEST_SINGLE_HIT, ";
    tmpStream << curTest.xtalIdx.getInt() << ", ";
    tmpStream << curTest.meV << ", ";
    tmpStream << curTest.xtalPos << ", ";

    if (curTest.trigMode == CalXtalId::BESTRANGE)
      tmpStream << "BESTRANGE, ";
    else tmpStream << "ALLRANGE, ";

    if (curTest.tupleOn) 
      tmpStream << "TUPLEON, " ;
    else tmpStream << "TUPLEOFF, ";

    if (curTest.gltOn)   
      tmpStream << "GLTON, " ;
    else tmpStream << "GLTOFF, ";

    if (curTest.noiseOn) 
      tmpStream << "NOISEON, " ;
    else tmpStream << "NOISEOFF, ";

    switch (curTest.getCalibSrc()) {
    case CurrentTest::FLIGHT_LAT:
      tmpStream << "FLIGHT_CALIB, ";
      break;
    case CurrentTest::PARTIAL_LAT:
      tmpStream << "8TOWER_CALIB, ";
      break;
    case CurrentTest::IDEAL:
      tmpStream << "IDEAL_CALIB, ";
      break;
    default:
      short src = curTest.getCalibSrc();
      assert(src >= 0 && src < CurrentTest::N_CALIB_SRC);
    }
    curTest.testDesc = tmpStream.str();
  }

  ////////////////////////////////////////////
  // SECTION 1: GENERATE MC INTEGRATING HIT //
  ////////////////////////////////////////////

  // used as output for xtalRecTuple
  CalTupleEntry calTupleEnt;

  Event::McIntegratingHit hit;
  fillMcHit(curTest.xtalIdx, 
            curTest.xtalPos, 
            curTest.meV,
            hit);

  //////////////////////////////
  //-- SECTION 2: XTAL DIGI --//
  //////////////////////////////
  // digi outputs
  Event::CalDigi calDigi(curTest.trigMode, curTest.xtalIdx.getCalXtalId()); 
  vector<const Event::McIntegratingHit*> hitVec;
  hitVec.push_back(&hit);
  Event::GltDigi glt;
  sc = curTest.getDigiTool()->calculate(hitVec,
                                        NULL,
                                        calDigi,
                                        curTest.lacBits,
                                        curTest.trigBits,
                                        (curTest.gltOn) ? &glt : NULL);
  if (sc.isFailure()) return sc;
  
  
  //-- VERIFY XTAL DIGI OUTPUT --//
  sc = verifyXtalDigi(calDigi, (curTest.gltOn) ? &glt : NULL);
  if (sc.isFailure()) return sc;

  
  // recon outputs
  Event::CalXtalRecData recData(CalXtalId::BESTRANGE, curTest.xtalIdx.getCalXtalId());
  CalArray<FaceNum, bool> belowThresh;
  bool xtalBelowThresh;
  CalArray<FaceNum, bool> saturated;

  ///////////////////////////////////
  //-- SECTION 3: RUN XTAL RECON --//
  ///////////////////////////////////
  sc = curTest.getRecTool()->calculate(calDigi,
                                       recData,
                                       belowThresh,
                                       xtalBelowThresh,
                                       saturated,
                                       curTest.tupleOn ? &calTupleEnt : 0);
  if (sc.isFailure()) return sc;

  // only process if ptr is non-null
  if (recData.getNReadouts() && !xtalBelowThresh) {
    Point testPoint;          // 3D point of test xtalPos
    pos2Point(curTest.xtalIdx, curTest.xtalPos, testPoint);
    curTest.posDiff = Vector(testPoint - recData.getPosition()).magnitude();
    curTest.recEne  = recData.getEnergy();
    float eneDiff = pct_diff(curTest.meV, curTest.recEne);

    //--------- REPORT RESULT --------------------------//
    ostringstream tmp;
    tmp << curTest.testDesc 
        << setfill(' ') << setw(7) << setprecision(2) << curTest.posDiff << ", "
        << setfill(' ') << setw(7) << fixed << setprecision(2) << eneDiff << "%";
    curTest.testDesc = tmp.str();
    msglog << MSG::INFO << curTest.testDesc << endreq;

    //-- skip this particular test if noise is on. b/c answers will vary too much
    if (!curTest.noiseOn) {
      //--------- TEST MARGINS ---------------------------//
      if (curTest.posDiff > m_singleHitPosMrgn) {
        msglog << MSG::WARNING << "TESTFAIL, BAD POS, "
               << curTest.posDiff << ", "
               << curTest.xtalPos << ", "
               << curTest.testDesc << endreq;
        //return StatusCode::FAILURE;
      }
      
      if (eneDiff > m_singleHitEneMrgn) {
        msglog << MSG::WARNING << "TESTFAIL, BAD ENE, "
               << eneDiff << ", "
               << curTest.meV << ", "
               << curTest.testDesc << endreq;
        //return StatusCode::FAILURE;
      }
    }

    //--------- OPTIONAL CAL TUPLE TEST -----------------//
    if (curTest.tupleOn) {
      // currently allways using 1st readout
      Event::CalDigi::CalXtalReadoutCol::const_iterator ro = 
        calDigi.getReadoutCol().begin();

      for (FaceNum face; face.isValid(); face++) {
        RngNum rng = ro->getRange(face);
        XtalRng xRng(face,rng);

        // tuple test1: adcPed
        float adcPedTpl = calTupleEnt.m_calXtalAdcPed[twr][lyr][col][face.getInt()];
        if (floor(adcPedTpl) != floor(curTest.adcPed[xRng])) {
          msglog << MSG::ERROR << "TESTFAIL, bad CALTUPLE.adcPed, "
                 << adcPedTpl << ", "
                 << curTest.adcPed[xRng] << ", "
                 << curTest.testDesc << endreq;
          return StatusCode::FAILURE;
        }
      }

      //-- these answers may vary too much w/ the noise on
      if (!curTest.noiseOn) {
        // tuple test2: faceSignal
        float fsigPos = calTupleEnt.m_calXtalFaceSignal[twr][lyr][col][POS_FACE];
        float fsigNeg = calTupleEnt.m_calXtalFaceSignal[twr][lyr][col][NEG_FACE];
        
        float meanSig = sqrt(fsigPos*fsigNeg);

        float fsigDiff = pct_diff(meanSig, curTest.meV);
        if ( fsigDiff > .05) {
          msglog << MSG::ERROR << "TESTFAIL, BAD CALTUPLE.faceSignal, "
                 << meanSig << ", "
                 << curTest.meV << ", "
                 << curTest.testDesc << endreq;
          return StatusCode::FAILURE;
        }
      }
    }
  }

  if (m_tuple) m_tuple->Fill();

  return StatusCode::SUCCESS;
}


StatusCode test_CalXtalResponse::lex8_to_lex1(FaceIdx faceIdx, float l8adc, float &l1adc, ICalCalibSvc &calCalibSvc) {

  float tmpDAC;
  StatusCode sc;
  MsgStream msglog(msgSvc(), name()); 

  FaceNum face = faceIdx.getFace();

  // use this point to generate LEX8/LEX1 ratio
  ValSig l8ULD;
  sc = calCalibSvc.getULDCI(RngIdx(faceIdx,LEX8), l8ULD);
  if (sc.isFailure()) return sc;

  // 1st convert to dac
  float x8tmp = curCalib.uldThresh[XtalRng(face,LEX8)];
  sc = calCalibSvc.evalDAC(RngIdx(faceIdx, LEX8),
                           x8tmp, tmpDAC);
  if (sc.isFailure()) return sc;
      
  // 2nd convert to LEX1 adc
  float x1tmp;
  sc = calCalibSvc.evalADC(RngIdx(faceIdx, LEX1),
                           tmpDAC, x1tmp);
  if (sc.isFailure()) return sc;
      
  float rat = x1tmp/x8tmp;

  l1adc = l8adc*rat;

  return StatusCode::SUCCESS;

}

StatusCode test_CalXtalResponse::hex8_to_hex1(FaceIdx faceIdx, float h8adc, float &h1adc, ICalCalibSvc &calCalibSvc) {

  float tmpDAC;
  StatusCode sc;
  MsgStream msglog(msgSvc(), name()); 

  FaceNum face = faceIdx.getFace();

  // use this point to generate HEX8/HEX1 ratio
  ValSig h8ULD;
  sc = calCalibSvc.getULDCI(RngIdx(faceIdx,HEX8), h8ULD);
  if (sc.isFailure()) return sc;

  // 1st convert to dac
  float x8tmp = curCalib.uldThresh[XtalRng(face,HEX8)];
  sc = calCalibSvc.evalDAC(RngIdx(faceIdx, HEX8),
                           x8tmp, tmpDAC);
  if (sc.isFailure()) return sc;
      
  // 2nd convert to HEX1 adc
  float x1tmp;
  sc = calCalibSvc.evalADC(RngIdx(faceIdx, HEX1),
                           tmpDAC, x1tmp);
  if (sc.isFailure()) return sc;
      
  float rat = x1tmp/x8tmp;

  h1adc = h8adc*rat;

  return StatusCode::SUCCESS;

}

StatusCode test_CalXtalResponse::verifyXtalDigi(const Event::CalDigi &calDigi,
                                                Event::GltDigi *glt) {
  StatusCode sc;
  MsgStream msglog(msgSvc(), name()); 

      
  //-- DIGI VERIFICTION: N READOUTS --//
  short nRO = calDigi.getReadoutCol().size();
  CalXtalId::CalTrigMode trigMode = calDigi.getMode();
  if ((trigMode == CalXtalId::BESTRANGE && nRO != 1) ||
      (trigMode == CalXtalId::ALLRANGE && nRO != 4)) {
    msglog << MSG::ERROR << "TESTFAIL, BAD N READOUTS, "
           << nRO << ", , "
           << curTest.testDesc << endreq;
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
    RngNum rng = ro->getRange(face);
    XtalRng xRng(face,rng);

    //////////////
    // LAC test //
    //////////////
    if (!trig_test_margin(curTest.fsignl[xRng], 
                          curCalib.lacMeV[face], 
                          curTest.lacBits[face],
                          curCalib.mevPerADC[xRng]) 
        ) {
      msglog << MSG::ERROR << "TESTFAIL, BAD LAC BIT, "
             << curTest.fsignl[xRng] << ", "
             << curCalib.lacMeV[face] << ", "
             << curTest.testDesc << endreq;
      return StatusCode::FAILURE;
    }
    

    //////////////
    // ULD test //
    //////////////
    //-- ULD test1: val should be < uld limit for current range
    // ULD test1: val should be < uld limit for selected range (w/in 1 ADC UNIT)
    if (!trig_test_margin(curTest.adcPed[xRng], curCalib.uldThresh[xRng], false, 1)) {
      msglog << MSG::ERROR << "TESTFAIL, BAD ULD, "
             << curTest.adcPed[xRng] << ", "
             << curCalib.uldThresh[xRng] << ", "
             << curTest.testDesc << endreq;
      return StatusCode::FAILURE;
    }
 
    // ULD test2: val should be > uld limit for next range down
    if (rng > LEX8) {
      if (!curTest.noiseOn) //noise can mess up this test b/c it differs per channel
        if (!trig_test_margin(curTest.fsignl[xRng], 
                              curCalib.uldMeV[XtalRng(face,rng.getInt()-1)], 
                              true, 
                              // xtra 1% margin for possible variation in faceSignal per channel
                              curCalib.mevPerADC[xRng] + .01*curTest.fsignl[xRng])) {
          msglog << MSG::ERROR << "TESTFAIL, BAD ULD, "
                 << curTest.fsignl[xRng] << ", "
                 << curCalib.uldMeV[XtalRng(face,rng.getInt()-1)] << ", "
                 << curTest.testDesc << endreq;
          return StatusCode::FAILURE;
        }
    }
  } // face loop

  sc = verifyXtalTrig(calDigi, curTest.trigBits, glt);
  if (sc.isFailure()) return sc;

  //-- CalTrigTool test --//
  CalArray<XtalDiode, bool> trigBitsTest;
  Event::GltDigi glt2;
  sc = curTest.getTrigTool()->calcXtalTrig(calDigi,
                                           trigBitsTest,
                                           (curTest.gltOn) ? &glt2 : NULL);
  if (sc.isFailure()) return sc;
  sc = verifyXtalTrig(calDigi, trigBitsTest, (curTest.gltOn) ? &glt2 : NULL);
  if (sc.isFailure()) return sc;

  //-- caltrigtool pass 2: this time w/ readout object & not calDigi object --//
  trigBitsTest.fill(false);
  Event::GltDigi glt3;
  sc = curTest.getTrigTool()->calcXtalTrig(curTest.xtalIdx,
                                           *ro,
                                           trigBitsTest,
                                           (curTest.gltOn) ? &glt3 : NULL);
  if (sc.isFailure()) return sc;
  sc = verifyXtalTrig(calDigi,trigBitsTest, (curTest.gltOn) ? &glt3 : NULL);
  if (sc.isFailure()) return sc;
                                           

  return StatusCode::SUCCESS;
}


StatusCode test_CalXtalResponse::verifyXtalTrig(const Event::CalDigi &calDigi,
                                                CalArray<XtalDiode, bool> &trigBits,
                                                Event::GltDigi *glt) {
  MsgStream msglog(msgSvc(), name()); 

  // currently allways using 1st readout
  Event::CalDigi::CalXtalReadoutCol::const_iterator ro = 
    calDigi.getReadoutCol().begin();


  for (FaceNum face; face.isValid(); face++) {
    FaceIdx faceIdx(curTest.xtalIdx, face);

    RngNum rng = ro->getRange(face);
    XtalRng xRng(face,rng);
          
                                    
    //////////////
    // FLE test.//
    //////////////

    //-- any test based on estimated faceSignal can suffer due to noise since
    //   noise varies from channel to channel
    if (!curTest.noiseOn)
      if (!trig_test_margin(curTest.fsignl[xRng], 
                            curCalib.trigMeV[XtalDiode(face,LRG_DIODE)], 
                            trigBits[XtalDiode(face,LRG_DIODE)], 
                            // xtra 1% margin for possible variation in faceSignal per channel
                            curCalib.mevPerADC[xRng]+.01*curTest.fsignl[xRng])) {
        msglog << MSG::ERROR << "TESTFAIL, BAD FLE BIT, "
               << curTest.fsignl[xRng] << ", "
               << curCalib.trigMeV[XtalDiode(face,LRG_DIODE)] << ", "
               << curTest.testDesc << endreq;
        return StatusCode::FAILURE;
      }
        
    //////////////
    // FHE test.//
    //////////////
    //-- any test based on estimated faceSignal can suffer due to noise since
    //   noise varies from channel to channel
    if (!curTest.noiseOn)
      if (!trig_test_margin(curTest.fsignl[xRng], 
                            curCalib.trigMeV[XtalDiode(face,SM_DIODE)], 
                            trigBits[XtalDiode(face,SM_DIODE)], 
                            // xtra 1% margin for possible variation in faceSignal per channel
                            curCalib.mevPerADC[xRng]+.01*curTest.fsignl[xRng])) {
        msglog << MSG::ERROR << "TESTFAIL, BAD FHE BIT, "
               << curTest.fsignl[xRng] << ", "
               << curCalib.trigMeV[XtalDiode(face,SM_DIODE)] << ", "
               << curTest.testDesc << endreq;
        //return StatusCode::FAILURE;
      }

    if (glt) {
      // test xtal FLE trigger
      if (trigBits[XtalDiode(face, LRG_DIODE)] !=
          glt->getCALLOtrigger(faceIdx.getCalXtalId())) {
        msglog << MSG::ERROR << "TESTFAIL, BAD GLT BIT, "
               << trigBits[XtalDiode(face, LRG_DIODE)] << ", "
               << glt->getCALLOtrigger(faceIdx.getCalXtalId()) << ", "
               << curTest.testDesc << endreq;
        return StatusCode::FAILURE;
      }
      
      // test xtal FHE trigger
      if (trigBits[XtalDiode(face, SM_DIODE)] !=
          glt->getCALHItrigger(faceIdx.getCalXtalId())) {
        msglog << MSG::ERROR << "TESTFAIL, BAD GLT BIT, "
               << trigBits[XtalDiode(face, SM_DIODE)] << ", "
               << glt->getCALHItrigger(faceIdx.getCalXtalId()) << ", "
               << curTest.testDesc << endreq;
        return StatusCode::FAILURE;
      }
    }
  }

  if (glt) {
    // test global FLE trigger
    if ((trigBits[XtalDiode(POS_FACE, LRG_DIODE)] ||
         trigBits[XtalDiode(NEG_FACE, LRG_DIODE)]) !=
        glt->getCALLOtrigger()) {
      msglog << MSG::ERROR << "TESTFAIL, BAD GLT BIT, "
             << (trigBits[XtalDiode(POS_FACE, LRG_DIODE)] ||
                 trigBits[XtalDiode(NEG_FACE, LRG_DIODE)]) << ", "
             << glt->getCALLOtrigger() << ", "
             << curTest.testDesc 
             << glt->getCALLOtrigger() << ", " << endreq;
      return StatusCode::FAILURE;
    }

    // test global FHE trigger
    if((trigBits[XtalDiode(POS_FACE, SM_DIODE)] ||
        trigBits[XtalDiode(NEG_FACE, SM_DIODE)]) != 
       glt->getCALHItrigger()) {
      msglog << MSG::ERROR << "TESTFAIL, BAD GLT BIT, "
             << (trigBits[XtalDiode(POS_FACE, SM_DIODE)] ||
                 trigBits[XtalDiode(NEG_FACE, SM_DIODE)]) << ", "
             << glt->getCALHItrigger() << ", "
             << curTest.testDesc << endreq;
      return StatusCode::FAILURE;
    }
  }

  return StatusCode::SUCCESS;
}


StatusCode test_CalXtalResponse::preprocXtalDigi(const Event::CalDigi &calDigi) {
  MsgStream msglog(msgSvc(), name());
  StatusCode sc;

  int nReadouts = (curTest.trigMode == CalXtalId::BESTRANGE) ? 1 : 4;


  for (int nRO = 0; nRO < nReadouts; nRO++) {
    for (FaceNum face; face.isValid(); face++) {
      float adc = calDigi.getAdc(nRO, face);
      RngNum rng = calDigi.getRange(nRO, face);
      RngIdx rngIdx(curTest.xtalIdx, face, rng);
      XtalRng xRng(face,rng);
      
      curTest.adcPed[xRng] = adc - curCalib.ped[xRng];


      sc = curTest.getCalCalibSvc()->evalFaceSignal(rngIdx, curTest.adcPed[xRng],
                                                    curTest.fsignl[xRng]);
      if (sc.isFailure()) return sc;
    }
  }


  return StatusCode::SUCCESS;
}


StatusCode test_CalXtalResponse::testCalCalibSvc() {
  StatusCode  sc;
  MsgStream msglog(msgSvc(), name()); 

  curCalib.clear();

  /////////////////////////////////////////////////
  // SECTION 0: GENERATE TEST DESCRIPTION STRING //
  /////////////////////////////////////////////////
  {
    ostringstream tmpStream;
    tmpStream << "TEST_CALCALIBSVC, ";
    tmpStream << curTest.xtalIdx.getInt() << ", ";

    switch (curTest.getCalibSrc()) {
    case CurrentTest::FLIGHT_LAT:
      tmpStream << "FLIGHT_CALIB, ";
      break;
    case CurrentTest::PARTIAL_LAT:
      tmpStream << "8TOWER_CALIB, ";
      break;
    case CurrentTest::IDEAL:
      tmpStream << "IDEAL_CALIB, ";
      break;
    default:
      short src = curTest.getCalibSrc();
      assert(src >= 0 && src < CurrentTest::N_CALIB_SRC);
    }
    curTest.testDesc = tmpStream.str();
  }


  //-- MPD --//
  ValSig mpdLrg, mpdSm;
  sc = curTest.getCalCalibSvc()->getMPD(curTest.xtalIdx,
                                        mpdLrg, mpdSm);
  if (sc.isFailure()) return sc;

  CalArray<DiodeNum, float> mpd;
  sc = curTest.getCalCalibSvc()->getMPD(curTest.xtalIdx, mpd);
  if (sc.isFailure()) return sc;

  if (mpd[LRG_DIODE] != mpdLrg.getVal()) {
    msglog << MSG::ERROR << "TESTFAIL, BAD MPD, "
           << curTest.testDesc << endreq;
    return StatusCode::FAILURE;
  }


  //-- PEDS --//
  sc = curTest.getCalCalibSvc()->getPed(curTest.xtalIdx,
                                        curCalib.ped, 
                                        curCalib.pedSig);
  if (sc.isFailure()) return sc;

  //-- ASYM --//
  const vector<CalibData::ValSig> *asymLrg;
  const vector<CalibData::ValSig> *asymSm;
  const vector<CalibData::ValSig> *asymNSPB;
  const vector<CalibData::ValSig> *asymPSNB;
  const vector<float>  *xVals;
  sc = curTest.getCalCalibSvc()->getAsym(curTest.xtalIdx,
                                         asymLrg,
                                         asymSm,
                                         asymNSPB,
                                         asymPSNB,
                                         xVals);
  if (sc.isFailure()) return sc;
  
  if (asymLrg->size()  <= 0 ||
      asymSm->size()   <= 0 ||
      asymNSPB->size() <= 0 ||
      asymPSNB->size() <= 0 ||
      xVals->size()     <= 0 ||
      xVals->size() != asymLrg->size() ||
      xVals->size() != asymSm->size() ||
      xVals->size() != asymPSNB->size() ||
      xVals->size() != asymNSPB->size()) {
    msglog << MSG::ERROR << "TESTFAIL, BAD ASYM, "
           << curTest.testDesc << endreq;
    return StatusCode::FAILURE;
  }

  for (unsigned i = 0; i < xVals->size(); i++) {
    float testAsym, testPos;

    sc = curTest.getCalCalibSvc()->evalPosLrg(curTest.xtalIdx, asymLrg->at(i).getVal(), testPos);
    if (sc.isFailure()) return sc;
    sc = curTest.getCalCalibSvc()->evalAsymLrg(curTest.xtalIdx, xVals->at(i), testAsym);
    if (sc.isFailure()) return sc;
    if (pct_diff(asymLrg->at(i).getVal(), testAsym) > .01 ||
        pct_diff(xVals->at(i), testPos) > .01) {
      msglog << MSG::ERROR << "TESTFAIL, BAD ASYM, "
             << curTest.testDesc << endreq;
      return StatusCode::FAILURE;
    }
    

    sc = curTest.getCalCalibSvc()->evalPosSm(curTest.xtalIdx, asymSm->at(i).getVal(), testPos);
    if (sc.isFailure()) return sc;
    sc = curTest.getCalCalibSvc()->evalAsymSm(curTest.xtalIdx, xVals->at(i), testAsym);
    if (sc.isFailure()) return sc;
    if (pct_diff(asymSm->at(i).getVal(), testAsym) > .01 ||
        pct_diff(xVals->at(i), testPos) > .01) {
      msglog << MSG::ERROR << "TESTFAIL, BAD ASYM, "
             << curTest.testDesc << endreq;
      return StatusCode::FAILURE;
    }

    sc = curTest.getCalCalibSvc()->evalPosNSPB(curTest.xtalIdx, asymNSPB->at(i).getVal(), testPos);
    if (sc.isFailure()) return sc;
    sc = curTest.getCalCalibSvc()->evalAsymNSPB(curTest.xtalIdx, xVals->at(i), testAsym);
    if (sc.isFailure()) return sc;
    if (pct_diff(asymNSPB->at(i).getVal(), testAsym) > .01 ||
        pct_diff(xVals->at(i), testPos) > .01) {
      msglog << MSG::ERROR << "TESTFAIL, BAD ASYM, "
             << curTest.testDesc << endreq;
      return StatusCode::FAILURE;
    }

    sc = curTest.getCalCalibSvc()->evalPosPSNB(curTest.xtalIdx, asymPSNB->at(i).getVal(), testPos);
    if (sc.isFailure()) return sc;
    sc = curTest.getCalCalibSvc()->evalAsymPSNB(curTest.xtalIdx, xVals->at(i), testAsym);
    if (sc.isFailure()) return sc;
    if (pct_diff(asymPSNB->at(i).getVal(), testAsym) > .01 ||
        pct_diff(xVals->at(i), testPos) > .01) {
      msglog << MSG::ERROR << "TESTFAIL, BAD ASYM, "
             << curTest.testDesc << endreq;
      return StatusCode::FAILURE;
    }
    
  } // xVals

    //-- THOLDCI --//
  sc = curTest.getCalCalibSvc()->getTholdCI(curTest.xtalIdx,
                                            curCalib.trigThresh,
                                            curCalib.lacThresh);
  if (sc.isFailure()) return sc;
  
  sc = curTest.getCalCalibSvc()->getULDCI(curTest.xtalIdx,
                                          curCalib.uldThresh);
  if (sc.isFailure()) return sc;
  
                                            

  //-- INTNONLIN && THOLDCI--//
  for (FaceNum face; face.isValid(); face++) {
    FaceIdx faceIdx(curTest.xtalIdx, face);
    ValSig fle, fhe, lac;
    sc = curTest.getCalCalibSvc()->getTholdCI(faceIdx, fle, fhe, lac);
    if (sc.isFailure()) return sc;
    
    if (fle.getVal() != curCalib.trigThresh[XtalDiode(face, LRG_DIODE)] ||
        fhe.getVal() != curCalib.trigThresh[XtalDiode(face, SM_DIODE)] ||
        lac.getVal() != curCalib.lacThresh[face]) {
      msglog << MSG::ERROR << "TESTFAIL, BAD THOLDCI, "
             << curTest.testDesc << endreq;
      return StatusCode::FAILURE;
    }

    for (RngNum rng; rng.isValid(); rng++) {
      XtalRng xRng(face,rng);
      RngIdx rngIdx(curTest.xtalIdx, xRng);
      
      ValSig uld;
      sc = curTest.getCalCalibSvc()->getULDCI(rngIdx, uld);
      if (sc.isFailure()) return sc;
      if (curCalib.uldThresh[xRng] != uld.getVal()) {
        msglog << MSG::ERROR << "TESTFAIL, BAD THOLDCI, "
               << curTest.testDesc << endreq;
        return StatusCode::FAILURE;
      }

      const vector<float> *adcs;
      const vector<float> *dacs;
      float err;

      sc = curTest.getCalCalibSvc()->getIntNonlin(rngIdx, adcs, dacs, err);
      if (sc.isFailure()) return sc;

      if (adcs->size() <=0 || 
          dacs->size() <= 0 || 
          adcs->size() > dacs->size()) {
        msglog << MSG::ERROR << "TESTFAIL, BAD INL, "
               << curTest.testDesc << endreq;
        return StatusCode::FAILURE;
      }

      for (unsigned i = 0; i < adcs->size(); i++) {
        float testDAC;
        sc = curTest.getCalCalibSvc()->evalDAC(rngIdx, adcs->at(i), testDAC);
        if (sc.isFailure()) return sc;

        if (pct_diff(testDAC, dacs->at(i)) > .01) {
          msglog << MSG::ERROR << "TESTFAIL, BAD INL, "
                 << curTest.testDesc << endreq;
          return StatusCode::FAILURE;
        }

        float testADC;
        sc = curTest.getCalCalibSvc()->evalADC(rngIdx, dacs->at(i), testADC);
        if (sc.isFailure()) return sc;

        if (pct_diff(testADC, adcs->at(i)) > .01) {
          msglog << MSG::ERROR << "TESTFAIL, BAD INL, "
                 << curTest.testDesc << endreq;
          return StatusCode::FAILURE;
        }

      } // inl points

      float ped, sig, cos;
      sc = curTest.getCalCalibSvc()->getPed(rngIdx, ped,sig,cos);
      if (sc.isFailure()) return sc;

      if (curCalib.ped[xRng] != ped ||
          curCalib.pedSig[xRng] != sig) {
        msglog << MSG::ERROR << "TESTFAIL, BAD PEDS, "
               << curTest.testDesc << endreq;
        return StatusCode::FAILURE;
      }
    }
  } // xRng

  return StatusCode::SUCCESS;
}


StatusCode test_CalXtalResponse::preprocCalCalib() {
  StatusCode  sc;
  MsgStream msglog(msgSvc(), name()); 
  
  for (FaceNum face; face.isValid(); face++) {
    FaceIdx faceIdx(curTest.xtalIdx, face);
    sc = curTest.getCalCalibSvc()->evalFaceSignal(RngIdx(faceIdx,LEX8),
                                                  curCalib.lacThresh[face],
                                                  curCalib.lacMeV[face]);
    if (sc.isFailure()) return sc;

    RngNum fleRng(LEX8);
    float tmpThresh = curCalib.trigThresh[XtalDiode(face,LRG_DIODE)];
    if (tmpThresh >= curCalib.uldThresh[XtalRng(face,LEX8)]) {
      fleRng = LEX1;
      sc = lex8_to_lex1(faceIdx, 
                        tmpThresh, 
                        tmpThresh, 
                        *curTest.getCalCalibSvc());
      if (sc.isFailure()) return sc;
    }
    sc = curTest.getCalCalibSvc()->evalFaceSignal(RngIdx(faceIdx,fleRng),
                                                  tmpThresh,
                                                  curCalib.trigMeV[XtalDiode(face, LRG_DIODE)]);
    if (sc.isFailure()) return sc;

    RngNum fheRng(HEX8);
    tmpThresh = curCalib.trigThresh[XtalDiode(face,SM_DIODE)];
    if (tmpThresh>= curCalib.uldThresh[XtalRng(face,HEX8)]) {
      fheRng = HEX1;
      sc = hex8_to_hex1(faceIdx, tmpThresh, tmpThresh, *curTest.getCalCalibSvc());
      if (sc.isFailure()) return sc;
    }
    sc = curTest.getCalCalibSvc()->evalFaceSignal(RngIdx(faceIdx,fheRng),
                                                  tmpThresh,
                                                  curCalib.trigMeV[XtalDiode(face, SM_DIODE)]);
    if (sc.isFailure()) return sc;
                                                  
    for (RngNum rng; rng.isValid(); rng++) {
      RngIdx rngIdx(faceIdx,rng);
      XtalRng xRng(face, rng);
      
      sc = curTest.getCalCalibSvc()->evalFaceSignal(rngIdx,
                                                    curCalib.uldThresh[xRng],
                                                    curCalib.uldMeV[xRng]);
      if (sc.isFailure()) return sc;
      
      // adcpermev - use last point in adc range
      curCalib.mevPerADC[xRng] = curCalib.uldMeV[xRng] / 
        curCalib.uldThresh[xRng];

      // calculate pedestal noise in MeV
      curCalib.pedNoise[xRng] = curCalib.pedSig[xRng] *
        curCalib.mevPerADC[xRng];
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
  curTest.xtalPos = 0;

  //-- HEX1 SATURATION --//
  // run meV up to 95% of saturation for lowest saturating face
  // should work fine.
  curTest.meV = min(curCalib.uldMeV[XtalRng(POS_FACE,HEX1)],
                    curCalib.uldMeV[XtalRng(NEG_FACE,HEX1)])*.95;
  sc = testSingleHit();
  if (sc.isFailure()) return sc;
  

  for (FaceNum face; face.isValid(); face++) {
    //-- LAC --//
    curTest.meV = curCalib.lacMeV[face]*.95;
    sc = testSingleHit();
    if (sc.isFailure()) return sc;

    curTest.meV = curCalib.lacMeV[face]*1.05;
    sc = testSingleHit();
    if (sc.isFailure()) return sc;

    for (DiodeNum diode; diode.isValid(); diode++) {
      XtalDiode xDiode(face,diode);
      //-- FLE/FHE --//
      curTest.meV = curCalib.trigMeV[xDiode]*.95;
      sc = testSingleHit();
      if (sc.isFailure()) return sc;

      curTest.meV = curCalib.trigMeV[xDiode]*1.05;
      sc = testSingleHit();
      if (sc.isFailure()) return sc;
      
    }

    //- HEX1 ULD is not a ULD, it's a saturation point.
    for (RngNum rng; rng < HEX1; rng++) {
      XtalRng xRng(face,rng);
      //-- ULD --//
      curTest.meV = curCalib.uldMeV[xRng]*.95;
      sc = testSingleHit();
      if (sc.isFailure()) return sc;

      curTest.meV = curCalib.uldMeV[xRng]*1.05;
      sc = testSingleHit();
      if (sc.isFailure()) return sc;
      
    }
  }

  return StatusCode::SUCCESS;
}

void test_CalXtalResponse::fillMcHit(XtalIdx xtalIdx,
                                     float xtalPos,
                                     float meV,
                                     Event::McIntegratingHit &hit) {

  //-- Create Volume Id (snagged from XtalRecTool::pos2Point()
  TwrNum twr = xtalIdx.getTwr();
  LyrNum lyr = xtalIdx.getLyr();
  ColNum col = xtalIdx.getCol();
  idents::VolumeIdentifier segmId;
  segmId.append(m_eLATTowers);
  segmId.append(twr.getRow());
  segmId.append(twr.getCol());
  segmId.append(m_eTowerCAL);
  segmId.append(lyr);
  segmId.append(lyr.getDir()); 
  segmId.append(col);
  segmId.append(m_eXtal);

  //-- calculate xtal segment #
  short nSeg = (short)floor(((xtalPos + m_CsILength/2)  // range 0->CsILen(mm)
                             / m_CsILength)             // range 0->1
                            * m_nCsISeg);               // range 0->nCsISeg
            
  //-- make sure we clip the edges to segments 0 && (nSeg - 1)
  nSeg = max<short>(nSeg,0);
  nSeg = min<short>(nSeg,m_nCsISeg-1);

  //  apply seg # to volId
  segmId.append(nSeg);

  //-- calc 1st moment (distance from segment ctr)
  float segCtr = (((float)nSeg+.5) // range 0.5 -> 11.5
                  / m_nCsISeg)     // range 0 -> 1
    * m_CsILength;   // range 0 -> csiLen
            
  float xDiff = xtalPos + m_CsILength/2 - segCtr;

  HepPoint3D mom1(xDiff,0,0);

  //-- Create McIntegratingHit
  hit.setVolumeID(segmId);
  // add single deposit w/ given mev & vector from center 
  // of segment.
  hit.addEnergyItem(meV, NULL, mom1);

}


StatusCode test_CalXtalResponse::testNoise() {
  StatusCode  sc;
  MsgStream msglog(msgSvc(), name()); 

  curTest.noiseOn = true;
  curTest.nHits = m_nNoiseHits;
  // test all ranges simultaneously
  curTest.trigMode = CalXtalId::ALLRANGE;
  curTest.xtalPos = 0;

  Event::McIntegratingHit hit;
  fillMcHit(curTest.xtalIdx, 
            curTest.xtalPos, 
            curTest.meV,
            hit);

  // for calculating stdev & mean
  CalArray<XtalRng, float> adcSum;
  CalArray<XtalRng, float> adcSumSq;

  adcSum.fill(0);
  adcSumSq.fill(0);
  
  for (int nHit = 0; nHit < curTest.nHits; nHit++) {
    sc = testSingleHit();
    if (sc.isFailure()) return sc;

    //-- for each range, sum up adc values
    for (XtalRng xRng; xRng.isValid(); xRng++) {
      adcSum[xRng] += curTest.adcPed[xRng];
      adcSumSq[xRng] += pow(curTest.adcPed[xRng],2);
    }
  }

  /////////////////////////////////////////////////
  // SECTION 0: GENERATE TEST DESCRIPTION STRING //
  /////////////////////////////////////////////////
  {
    ostringstream tmpStream;
    tmpStream << "TEST_SINGLE_HIT, ";
    tmpStream << curTest.xtalIdx.getInt() << ", ";
    tmpStream << curTest.meV << ", ";

    switch (curTest.getCalibSrc()) {
    case CurrentTest::FLIGHT_LAT:
      tmpStream << "FLIGHT_CALIB, ";
      break;
    case CurrentTest::PARTIAL_LAT:
      tmpStream << "8TOWER_CALIB, ";
      break;
    case CurrentTest::IDEAL:
      tmpStream << "IDEAL_CALIB, ";
      break;
    default:
      short src = curTest.getCalibSrc();
      assert(src >= 0 && src < CurrentTest::N_CALIB_SRC);
    }
    curTest.testDesc = tmpStream.str();
  }


  // calc & check 
  for (XtalRng xRng; xRng.isValid(); xRng++) {
    // do not look at ranges for which the ene level is too high
    if (curTest.meV > curCalib.uldMeV[xRng]*.90) continue;
    // do not look at ranges for which the ene level is too low
    if (curTest.meV < curCalib.pedNoise[xRng]) continue;

    float adcMean = adcSum[xRng]/curTest.nHits;
    float meanSq  = adcSumSq[xRng]/curTest.nHits;
    float adcSig  = sqrt(meanSq - pow(adcMean,2));

    if (pct_diff(adcSig, curCalib.pedSig[xRng]) > .1) {
      msglog << MSG::WARNING << "TESTFAIL, BAD NOISE SIGMA, "
             << adcSig << ", " 
             << curCalib.pedSig[xRng] << ", "
             << curTest.testDesc << endreq;
      //return StatusCode::FAILURE;
    }

    float fsignlMean;
    sc = curTest.getCalCalibSvc()->evalFaceSignal(RngIdx(curTest.xtalIdx,xRng), adcMean, fsignlMean);
    if (pct_diff(fsignlMean, curTest.meV) > .05) {
      msglog << MSG::WARNING << "TESTFAIL, BAD NOISE MEAN, "
             << fsignlMean << ", "
             << curTest.meV << ", "
             << curTest.testDesc 
             << xRng.getFace().getInt() << ", "
             << xRng.getRng().getInt()  << ", " << endreq;
      //return StatusCode::FAILURE;
    }
  }

  return StatusCode::SUCCESS;
}


StatusCode test_CalXtalResponse::testMultiHit() {
  StatusCode  sc;
  MsgStream msglog(msgSvc(), name()); 

  curTest.nHits   = 10;
  curTest.gltOn   = false;
  curTest.noiseOn = false;
  curTest.tupleOn = false;

  int nPosPts = 5;



  //////////////////////
  // TEST A: IN-PLACE //
  //////////////////////
  // test several hits in same part of xtal.
  // perform this test at 5 places along xtal

  for (int nPos = 0; nPos < nPosPts; nPos++) {
    // longitudinal distance from ctr of xtal in mm
    curTest.xtalPos = m_CsILength*nPos/(nPosPts-1) - m_CsILength/2;

    // to hold actual mc hits
    vector<Event::McIntegratingHit> hitVec;
    
    // create individual hits
    for (int nHit = 0; nHit < curTest.nHits; nHit++) {
      Event::McIntegratingHit hit;
      fillMcHit(curTest.xtalIdx, 
                curTest.xtalPos, 
                curTest.meV/curTest.nHits,
                hit);
      hitVec.push_back(hit);
    }

    /////////////////////////////////////////////////
    // SECTION 0: GENERATE TEST DESCRIPTION STRING //
    /////////////////////////////////////////////////
    {
      ostringstream tmpStream;
      tmpStream << "TEST_MULTIHIT_SINGLEPOS, ";
      tmpStream << curTest.xtalIdx.getInt() << ", ";
      tmpStream << curTest.meV << ", ";
      tmpStream << curTest.xtalPos << ", ";

      switch (curTest.getCalibSrc()) {
      case CurrentTest::FLIGHT_LAT:
        tmpStream << "FLIGHT_CALIB, ";
        break;
      case CurrentTest::PARTIAL_LAT:
        tmpStream << "8TOWER_CALIB, ";
        break;
      case CurrentTest::IDEAL:
        tmpStream << "IDEAL_CALIB, ";
        break;
      default:
        short src = curTest.getCalibSrc();
        assert(src >= 0 && src < CurrentTest::N_CALIB_SRC);
      }
      curTest.testDesc = tmpStream.str();
    }

    //-- XTAL DIGI --//
    Event::CalDigi calDigi(curTest.trigMode, curTest.xtalIdx.getCalXtalId()); 
    // list of ptrs is needed for function call
    vector<const Event::McIntegratingHit*> hitVecPtr;
    for (unsigned i = 0; i < hitVec.size(); i++)
      hitVecPtr.push_back(&(hitVec[i]));
    sc = curTest.getDigiTool()->calculate(hitVecPtr,
                                          NULL,
                                          calDigi,
                                          curTest.lacBits,
                                          curTest.trigBits,
                                          NULL);
    if (sc.isFailure()) return sc;

    //-- VERIFY XTAL DIGI OUTPUT --//
    sc = verifyXtalDigi(calDigi, NULL);
    if (sc.isFailure()) return sc;

    // recon outputs
    Event::CalXtalRecData recData(CalXtalId::BESTRANGE, curTest.xtalIdx.getCalXtalId());
    CalArray<FaceNum, bool> belowThresh;
    bool xtalBelowThresh;
    CalArray<FaceNum, bool> saturated;

    ///////////////////////////////////
    //-- SECTION 3: RUN XTAL RECON --//
    ///////////////////////////////////
    sc = curTest.getRecTool()->calculate(calDigi,
                                         recData,
                                         belowThresh,
                                         xtalBelowThresh,
                                         saturated,
                                         0);
    if (sc.isFailure()) return sc;

    // only process if ptr is non-null
    if (recData.getNReadouts() && !xtalBelowThresh) {
      Point testPoint;          // 3D point of test xtalPos
      pos2Point(curTest.xtalIdx, curTest.xtalPos, testPoint);
      curTest.posDiff = Vector(testPoint - recData.getPosition()).magnitude();
      curTest.recEne  = recData.getEnergy();
      float eneDiff = pct_diff(curTest.meV, curTest.recEne);
        
      //--------- REPORT RESULT --------------------------//
      ostringstream tmp;
      tmp << curTest.testDesc 
          << setfill(' ') << setw(7) << setprecision(2) << curTest.posDiff << ", "
          << setfill(' ') << setw(7) << fixed << setprecision(2) << eneDiff << "%";
      curTest.testDesc = tmp.str();
      msglog << MSG::INFO << curTest.testDesc << endreq;
        
      //--------- TEST MARGINS ---------------------------//
      if (curTest.posDiff > m_singleHitPosMrgn) {
        msglog << MSG::ERROR << "TESTFAIL, BAD POS, "
               << curTest.posDiff << ", "
               << curTest.xtalPos << ", "
               << curTest.testDesc << endreq;
        //return StatusCode::FAILURE;
      }
      
      if (eneDiff > m_singleHitEneMrgn) {
        msglog << MSG::ERROR << "TESTFAIL, BAD ENE, "
               << eneDiff << ", "
               << curTest.meV << ", "
               << curTest.testDesc << endreq;
        //return StatusCode::FAILURE;
      }
    }
  }

  ////////////////////////
  // TEST B: ALONG XTAL //
  ////////////////////////
  // test several hits spread out along xtal
  curTest.xtalPos = 0; // average position will allways be at ctr 


  // to hold actual mc hits
  vector<Event::McIntegratingHit> hitVec;

  // create individual hits
  for (int nHit = 0; nHit < curTest.nHits; nHit++) {
    float xtalPos = m_CsILength*nHit/(curTest.nHits-1) - m_CsILength/2;
    Event::McIntegratingHit hit;
    fillMcHit(curTest.xtalIdx, 
              xtalPos, 
              curTest.meV/curTest.nHits,
              hit);
    hitVec.push_back(hit);    
  }

  /////////////////////////////////////////////////
  // SECTION 0: GENERATE TEST DESCRIPTION STRING //
  /////////////////////////////////////////////////
  {
    ostringstream tmpStream;
    tmpStream << "TEST_MULTIHIT_MULTIPOS, ";
    tmpStream << curTest.xtalIdx.getInt() << ", ";
    tmpStream << curTest.meV << ", ";

    switch (curTest.getCalibSrc()) {
    case CurrentTest::FLIGHT_LAT:
      tmpStream << "FLIGHT_CALIB, ";
      break;
    case CurrentTest::PARTIAL_LAT:
      tmpStream << "8TOWER_CALIB, ";
      break;
    case CurrentTest::IDEAL:
      tmpStream << "IDEAL_CALIB, ";
      break;
    default:
      short src = curTest.getCalibSrc();
      assert(src >= 0 && src < CurrentTest::N_CALIB_SRC);
    }
    curTest.testDesc = tmpStream.str();
  }

  //-- XTAL DIGI --//
  Event::CalDigi calDigi(curTest.trigMode, curTest.xtalIdx.getCalXtalId()); 
  
  // list of ptrs is needed for function call
  vector<const Event::McIntegratingHit*> hitVecPtr;
  for (unsigned i = 0; i < hitVec.size(); i++)
    hitVecPtr.push_back(&(hitVec[i]));
  sc = curTest.getDigiTool()->calculate(hitVecPtr,
                                        NULL,
                                        calDigi,
                                        curTest.lacBits,
                                        curTest.trigBits,
                                        NULL);
  if (sc.isFailure()) return sc;

  //-- VERIFY XTAL DIGI OUTPUT --//
  sc = verifyXtalDigi(calDigi, NULL);
  if (sc.isFailure()) return sc;

  // recon outputs
  Event::CalXtalRecData recData(CalXtalId::BESTRANGE, curTest.xtalIdx.getCalXtalId());
  CalArray<FaceNum, bool> belowThresh;
  bool xtalBelowThresh;
  CalArray<FaceNum, bool> saturated;

  ///////////////////////////////////
  //-- SECTION 3: RUN XTAL RECON --//
  ///////////////////////////////////
  sc = curTest.getRecTool()->calculate(calDigi,
                                       recData,
                                       belowThresh,
                                       xtalBelowThresh,
                                       saturated,
                                       0);
  if (sc.isFailure()) return sc;

  // only process if ptr is non-null
  if (recData.getNReadouts() && !xtalBelowThresh) {
    Point testPoint;          // 3D point of test xtalPos
    pos2Point(curTest.xtalIdx, curTest.xtalPos, testPoint);
    curTest.posDiff = Vector(testPoint - recData.getPosition()).magnitude();
    curTest.recEne  = recData.getEnergy();
    float eneDiff = pct_diff(curTest.meV, curTest.recEne);
        
    //--------- REPORT RESULT --------------------------//
    ostringstream tmp;
    tmp << curTest.testDesc 
        << setfill(' ') << setw(7) << setprecision(2) << curTest.posDiff << ", "
        << setfill(' ') << setw(7) << fixed << setprecision(2) << eneDiff << "%";
    curTest.testDesc = tmp.str();
    msglog << MSG::INFO << curTest.testDesc << endreq;
        
    //--------- TEST MARGINS ---------------------------//
    if (curTest.posDiff > m_singleHitPosMrgn) {
      msglog << MSG::ERROR << "TESTFAIL, BAD POS, "
             << curTest.posDiff << ", "
             << curTest.xtalPos << ", "
             << curTest.testDesc << endreq;
      //return StatusCode::FAILURE;
    }
      
    if (eneDiff > m_singleHitEneMrgn) {
      msglog << MSG::ERROR << "TESTFAIL, BAD ENE, "
             << eneDiff << ", "
             << curTest.meV << ", "
             << curTest.testDesc << endreq;
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

  /// inject 5% more energy than highest HEX1 saturation point for both faces.
  curTest.meV = 1.05*max(curCalib.uldMeV[XtalRng(POS_FACE, HEX1)],
                         curCalib.uldMeV[XtalRng(NEG_FACE, HEX1)]);
  // face signals are measured at center of xtal, so that's what we have to use.
  curTest.xtalPos = 0; 
  curTest.nHits = 1;
  curTest.gltOn = false;
  curTest.noiseOn = false;
  curTest.tupleOn = false;
  

  /////////////////////////////////////////////////
  // SECTION 0: GENERATE TEST DESCRIPTION STRING //
  /////////////////////////////////////////////////
  {
    ostringstream tmpStream;
    tmpStream << "TEST_HEX1_SATURATION, ";
    tmpStream << curTest.xtalIdx.getInt() << ", ";
    tmpStream << curTest.meV << ", ";

    switch (curTest.getCalibSrc()) {
    case CurrentTest::FLIGHT_LAT:
      tmpStream << "FLIGHT_CALIB, ";
      break;
    case CurrentTest::PARTIAL_LAT:
      tmpStream << "8TOWER_CALIB, ";
      break;
    case CurrentTest::IDEAL:
      tmpStream << "IDEAL_CALIB, ";
      break;
    default:
      short src = curTest.getCalibSrc();
      assert(src >= 0 && src < CurrentTest::N_CALIB_SRC);
    }
    curTest.testDesc = tmpStream.str();
  }

  Event::McIntegratingHit hit;
  fillMcHit(curTest.xtalIdx, 
            curTest.xtalPos, 
            curTest.meV,
            hit);


  //////////////////////////////
  //-- SECTION 2: XTAL DIGI --//
  //////////////////////////////
  // digi outputs
  Event::CalDigi calDigi(curTest.trigMode, curTest.xtalIdx.getCalXtalId()); 
  vector<const Event::McIntegratingHit*> hitVec;
  hitVec.push_back(&hit);
  sc = curTest.getDigiTool()->calculate(hitVec,
                                        NULL,
                                        calDigi,
                                        curTest.lacBits,
                                        curTest.trigBits,
                                        NULL);
  if (sc.isFailure()) return sc;
  
  
  //-- VERIFY XTAL DIGI OUTPUT --//
  sc = verifyXtalDigi(calDigi, NULL);
  if (sc.isFailure()) return sc;

  return StatusCode::SUCCESS;
}
