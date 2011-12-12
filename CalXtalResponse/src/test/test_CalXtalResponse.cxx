// $Header$

/** @file
    @author Z.Fewtrell
    
    @brief  CalXtalRepsonse unit test app - attempt to test as much as possible
    
*/

// LOCAL
#include "TestCfg.h"
#include "test_CalCalibSvc.h"
#include "test_CalSignalTool.h"
#include "test_CalDigi.h"
#include "test_CalXtalRecAlg.h"
#include "test_CalTrigTool.h"
#include "test_CalDiagnosticTool.h"
#include "test_PrecalcCalibTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"
#include "CalXtalResponse/ICalSignalTool.h"
#include "CalXtalResponse/ICalTrigTool.h"
#include "CalXtalResponse/ICalDiagnosticTool.h"
#include "../CalCalib/IPrecalcCalibTool.h"
#include "test_util.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "facilities/Util.h" // expandEnvVar()
#include "CalUtil/CalGeom.h" // findActiveTowers()
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Digi/CalDigi.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"

// EXTLIB
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"

// STD
#include <map>
#include <memory>


class ICalCalibSvc;
class IPrecalcCalibTool;
class ICalSignalTool;  
class ICalTrigTool;
class ICalDiagnosticTool;

/** unit test for CalXtalResponse package
    
    Structure:
    Individual unit test classes created for most major classes in
    CalXtalResponse.

    Algorithm:
    - A sequence of events goes through full MC->Digi->Recon process.
    - Each event has randomly selected energy deposits and jobOptions.
    - After each step, the appropriate sub-test is called for that
    particular tool.

    Variation of input parameters:
    - energy range from .1 MeV -> 100 GeV (logarithmic distribution to
    allow for equal # of hits in each energy range and more hits near
    'interesting' points like LAC/FLE/FHE/ULD thresholds
    - randomly selected deposit position along crystal length
    - multiple deposits per crystal
    - zero deposits per crystal
    - randomly selected crystals
    - multiple crystals
    - zeroSuppression on / off
    - BESTRANGE / ALLRANGE
    - ideal calibrations , real calibrations
    - partial LAT configurations (code still works with missing tower
    modules)

    Component Tests:
    - CalCalibSvc
    - PrecalcCalibSvc
    - CalSignalTool (incl XtalSignalTool)
    - CalDigiAlg (incl XtalDigiTool)
    - CalTrigTool (from MC)
    - CalTrigTool (from Digi)
    - CalDiagnosticTool
    - CalXtalRecAlg (incl (XtalRecTool)

    @ todo
    - test andrey's forced first range code (EZ)
    - direct diode deposits 
    - CalTupleAlg - not sure how to arrange this without creating
    unwanted output file.
    - noise simulation test
    - missing TDS collections (not just empty)
    - using TrigConfigSvc for cal readout mode selection
*/
class test_CalXtalResponse : 
  public Algorithm {
public:
  test_CalXtalResponse(const std::string& name, ISvcLocator* pSvcLocator);

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize(); 

private:
  StatusCode retrieveConstants(IGlastDetSvc &detSvc);

  /// apply current configuration to child tools & algorithms.
  StatusCode applyNewCfg(const CalXtalResponse::TestCfg &testCfg);

  /// create MC data consistent with current configuration.
  StatusCode fillMC(const CalXtalResponse::TestCfg &testCfg,
                    const CalXtalResponse::TwrSet &twrSet);

  /// verify CalCalibSvc
  StatusCode verifyCalCalibSvc(const CalXtalResponse::TestCfg &testCfg);

  /// fill mcIntegrating hit w/ given xtal position & energy
  /// \param xtalPos in mm along xtal length. from 0 - 326
  void fillMcHit(const CalUtil::XtalIdx xtalIdx,
                 const float xtalPos,
                 const float meV,
                 Event::McIntegratingHit &hit);

  /// register McIntegratingHit col with event service.  remove &
  /// delete any existing hit col
  StatusCode registerMcHitCol(Event::McIntegratingHitCol *hitCol);

  /// remove MCHits from TDS.
  StatusCode clearMC();

  /// need to be able to set properties on the fly for this alg
  /// must make it a child alg of current alg
  Algorithm *m_calDigiAlg;

  /// need to be able to set properties on the fly for this alg
  /// must make it a child alg of current alg
  Algorithm *m_calXtalRecAlg;


  ICalCalibSvc *m_calCalibSvc;

  /// CalCalibSvc w/ 'ideal mode' option enabled.
  ICalCalibSvc *m_calCalibSvcIdeal;

  ICalSignalTool *m_calSignalTool;

  /// used to calculate Cal Trig response from MC (default)
  ICalTrigTool *m_calTrigTool;

  /// used to calculate Cal Trig response from Digi
  ICalTrigTool *m_calTrigToolDigi;

  ICalDiagnosticTool *m_calDiagnosticTool;

  IPrecalcCalibTool *m_precalcCalibTool;

  /// list of active tower bays, populated at run time
  CalXtalResponse::TwrSet m_twrSet;

  /// list of all towers (used by ideal mode calib)
  static const CalXtalResponse::TwrSet m_fullTwrSet;

  /// current test case spec
  CalXtalResponse::TestCfg m_testCfg;

  /// read in calibrations from TXT file to validate against system
  /// provided values
  test_CalCalibSvc::TestCalibSet m_testCalibSet;

  /// read in ideal calibrations from TXT file to validate against system
  /// provided values
  test_CalCalibSvc::TestCalibSet m_testCalibSetIdeal;

  test_CalCalibSvc m_testCalCalibSvc;

  test_PrecalcCalibTool m_testPrecalcCalibTool;

  test_CalSignalTool m_testCalSignalTool;

  test_CalDigi m_testCalDigi;

  test_CalTrigTool m_testCalTrigTool;

  test_CalDiagnosticTool m_testCalDiagnosticTool;

  test_CalXtalRecAlg m_testCalXtalRecAlg;


  /// length of cal crystal in mm
  float m_csiLength;

  /// the value of fLATObjects field, defining LAT towers 
  int m_eLATTowers; 
  
  /// the value of fTowerObject field, defining calorimeter module 
  int m_eTowerCAL;
  /// the value of fCellCmp field defining CsI crystal
  int m_eXtal;      
  /// number of geometric segments per Xtal
  int m_nCsISeg; 

  /// path to txt table of calibration test values
  StringProperty m_pedTXTPath;
  /// path to txt table of calibration test values
  StringProperty m_cidac2adcTXTPath;
  /// path to txt table of calibration test values
  StringProperty m_mpdTXTPath;
  /// path to txt table of calibration test values
  StringProperty m_asymTXTPath;
  /// path to txt table of calibration test values
  StringProperty m_tholdCITXTPath;

  /// path to txt table of calibration test values
  StringProperty m_pedTXTPathIdeal;
  /// path to txt table of calibration test values
  StringProperty m_cidac2adcTXTPathIdeal;
  /// path to txt table of calibration test values
  StringProperty m_mpdTXTPathIdeal;
  /// path to txt table of calibration test values
  StringProperty m_asymTXTPathIdeal;
  /// path to txt table of calibration test values
  StringProperty m_tholdCITXTPathIdeal;

  /// keep count of events
  unsigned m_eventId;
  
}; // class test_CalXtalResponse

//static const AlgFactory<test_CalXtalResponse>  Factory;
//const IAlgFactory& test_CalXtalResponseFactory = Factory;
DECLARE_ALGORITHM_FACTORY(test_CalXtalResponse);

namespace {
  static const CalUtil::TwrNum _fullTwrSet[] = {0, 1, 2, 3, 
                                                4, 5, 6, 7, 
                                                8, 9, 10,11,
                                                12,13,14,15};
  
  /// expand environment variables in string without modifying
  /// original string
  std::string expandEnvVarConst(const std::string &input) {
    std::string copy(input);
    
    facilities::Util::expandEnvVar(&copy);

    return copy;
  }
  
}

using namespace CalXtalResponse;
using namespace std;
using namespace idents;
using namespace CalUtil;

const CalXtalResponse::TwrSet test_CalXtalResponse::m_fullTwrSet(_fullTwrSet,
                                                                 _fullTwrSet+16);

test_CalXtalResponse::test_CalXtalResponse(const string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator),
  m_calDigiAlg(0),
  m_calXtalRecAlg(0),
  m_calCalibSvc(0),
  m_calCalibSvcIdeal(0),
  m_calSignalTool(0),
  m_calTrigTool(0),
  m_calTrigToolDigi(0),
  m_calDiagnosticTool(0),
  m_precalcCalibTool(0),
  m_csiLength(0),
  m_eventId(0)
{
  declareProperty("pedTXTPath", m_pedTXTPath="");
  declareProperty("cidac2adcTXTPath", m_cidac2adcTXTPath="");
  declareProperty("mpdTXTPath", m_mpdTXTPath="");
  declareProperty("asymTXTPath", m_asymTXTPath="");
  declareProperty("tholdCITXTPath", m_tholdCITXTPath="");
 
  declareProperty("pedTXTPathIdeal", m_pedTXTPathIdeal="");
  declareProperty("cidac2adcTXTPathIdeal", m_cidac2adcTXTPathIdeal="");
  declareProperty("mpdTXTPathIdeal", m_mpdTXTPathIdeal="");
  declareProperty("asymTXTPathIdeal", m_asymTXTPathIdeal="");
  declareProperty("tholdCITXTPathIdeal", m_tholdCITXTPathIdeal="");

}

StatusCode test_CalXtalResponse::initialize() {
  StatusCode sc;
  MsgStream msglog(msgSvc(), name());   

  sc = createSubAlgorithm("CalDigiAlg", "CalDigiAlg", m_calDigiAlg);
  if (sc.isFailure())
    return sc;

  sc = createSubAlgorithm("CalXtalRecAlg", "CalXtalRecAlg", m_calXtalRecAlg);
  if (sc.isFailure())
    return sc;


  // obtain CalCalibSvc
  //sc = service("CalCalibSvc", "CalCalibSvc", m_calCalibSvc);
  sc = service("CalCalibSvc", m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get CalCalibSvc" << endreq;
    return sc;
  }

  //sc = service("CalCalibSvc", "CalCalibSvcIdeal", m_calCalibSvcIdeal);
  sc = service("CalCalibSvcIdeal", m_calCalibSvcIdeal);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get CalCalibSvcIdeal" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("CalSignalTool",
                               "CalSignalTool",
                               m_calSignalTool,
                               0); // intended to be shared
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create CalSignalTool" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("CalTrigTool",
                               "CalTrigTool",
                               m_calTrigTool,
                               0); // intended to be shared
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create CalTrigTool" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("CalTrigTool",
                               "CalTrigToolDigi",
                               m_calTrigToolDigi,
                               this); // not intended to be shared
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create CalTrigTool" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("CalDiagnosticTool",
                               "CalDiagnosticTool",
                               m_calDiagnosticTool,
                               0); // intended to be shared
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create CalDiagnosticTool" << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("PrecalcCalibTool",
                               "PrecalcCalibTool",
                               m_precalcCalibTool,
                               0); // intended to be shared
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create PrecalcCalibTool" << endreq;
    return sc;
  }

  // try to find the GlastDetSvc service
  IGlastDetSvc *detSvc = 0;
  sc = service("GlastDetSvc", detSvc);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't get GlastDetSvc " << endreq;
    return sc;
  }

  //-- find out which tems are installed.
  CalUtil::TwrList twrList = CalUtil::findActiveTowers(*detSvc);
  m_twrSet.insert(twrList.begin(), twrList.end());
  

  // initialize test_CalCalibSvc
  if (m_testCalCalibSvc.initialize(msgSvc()).isFailure())
    return StatusCode::FAILURE;

  if (m_testPrecalcCalibTool.initialize(msgSvc()).isFailure())
    return StatusCode::FAILURE;

  if (m_testCalSignalTool.initialize(msgSvc()).isFailure())
    return StatusCode::FAILURE;

  if (m_testCalDigi.initialize(msgSvc()).isFailure())
    return StatusCode::FAILURE;

  if (m_testCalTrigTool.initialize(msgSvc()).isFailure())
    return StatusCode::FAILURE;

  if (m_testCalDiagnosticTool.initialize(msgSvc()).isFailure())
    return StatusCode::FAILURE;

  if (m_testCalXtalRecAlg.initialize(msgSvc(), detSvc).isFailure())
    return StatusCode::FAILURE;

  sc = retrieveConstants(*detSvc);
  if (sc.isFailure())
    return sc;

  /// load calibration 'answers' from txt file
  m_testCalibSet.m_calPed.readTXT(expandEnvVarConst(m_pedTXTPath));
  m_testCalibSet.m_cidac2adc.readTXT(expandEnvVarConst(m_cidac2adcTXTPath));
  m_testCalibSet.m_calAsym.readTXT(expandEnvVarConst(m_asymTXTPath));
  m_testCalibSet.m_calMPD.readTXT(expandEnvVarConst(m_mpdTXTPath));
  m_testCalibSet.m_calTholdCI.readTXT(expandEnvVarConst(m_tholdCITXTPath));

  m_testCalibSetIdeal.m_calPed.readTXT(expandEnvVarConst(m_pedTXTPathIdeal));
  m_testCalibSetIdeal.m_cidac2adc.readTXT(expandEnvVarConst(m_cidac2adcTXTPathIdeal));
  m_testCalibSetIdeal.m_calAsym.readTXT(expandEnvVarConst(m_asymTXTPathIdeal));
  m_testCalibSetIdeal.m_calMPD.readTXT(expandEnvVarConst(m_mpdTXTPathIdeal));
  m_testCalibSetIdeal.m_calTholdCI.readTXT(expandEnvVarConst(m_tholdCITXTPathIdeal));

  
  return StatusCode::SUCCESS;
}

/// retrieve instrument constants from database
StatusCode test_CalXtalResponse::retrieveConstants(IGlastDetSvc &detSvc) {
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
    if(!detSvc.getNumericConstByName((*iter).second, &tmp)) {
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::ERROR << " constant " <<(*iter).second
             << " not defined" << endreq;
      return StatusCode::FAILURE;
    } else *((*iter).first)= int(tmp); // store retrieved value 
  }

  // DOUBLE CONSTANTS
  StatusCode sc;
  sc = detSvc.getNumericConstByName("CsILength", &tmp);
  m_csiLength = tmp;
  if (sc.isFailure()) {
    MsgStream msglog(msgSvc(), name());
    msglog << MSG::ERROR << " constant CsILength not defined" << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}

/**
   Alogrithm:
   For each new event:
   - generate random test case configuration
   - apply new cfg to any Gaudi child algorithms & tools (alter their
   jobOptions)
   - verify CalCalibSvc
   - verify precalcCalibTool
   - generate simulated MCInteratingHit deposit records from TestCfg
   - verify CalSignalTool (CalSignalTool should  process all the new
   MC data on the first query of the event.
   - execute CalDigiAlg
   - verify CalDigiCol (output of CalDigiAlg)
   - verify CalTrigTool (it should calculate from MC data since it is
   present in TDS
   - remove MC data from TDS (force CalTrigTool to use CalDigis)
   - verify CalTrigTool again (this time with CalDigi data instead of
   MC)
   - verify CalDiagnosticTool
   - run CalXtalRecAlg
   - verify CalXtalRecCOl (output of CalXtalRecAlg)
   
*/
StatusCode test_CalXtalResponse::execute() {
  StatusCode sc;
  MsgStream msglog(msgSvc(), name());  

  m_eventId++;

  /// generate new random test event.
  m_testCfg.randomFill();
  msglog << MSG::INFO << "event=" << m_eventId << " " << m_testCfg << endreq;

  /// apply new random configuration to Gaudi child algorithms & tools
  sc = applyNewCfg(m_testCfg);
  if (sc.isFailure()) 
    return sc;

  /// check calibrations
  msglog << MSG::DEBUG << "verifyCalCalibSvc" << endreq;
  sc = verifyCalCalibSvc(m_testCfg);
  if (sc.isFailure())
    return sc;

  /// check dreived calibrations
  msglog << MSG::DEBUG << "verify precalcCalibTool" << endreq;
  if (m_testPrecalcCalibTool.verify(*m_precalcCalibTool, 
                                    *m_calCalibSvc,
                                    m_testCfg, 
                                    m_twrSet).isFailure())
    return StatusCode::FAILURE;

  /// fill MC
  sc = fillMC(m_testCfg, m_twrSet);
  if (sc.isFailure())
    return sc;

  // check CalSignalTool
  msglog << MSG::DEBUG << "test CalSignalTool" << endreq;
  sc = m_testCalSignalTool.verify(*m_calSignalTool,
                                  *m_calCalibSvc,
                                  m_testCfg,
                                  m_twrSet);
  if (sc.isFailure())
    return sc;
                                 

  // execute CalDigiAlg
  sc = m_calDigiAlg->execute();
  if (sc.isFailure())
    return sc;

  // retrieve caldigicol
  // get a pointer to the input TDS data collection
  SmartDataPtr<Event::CalDigiCol> calDigiCol = SmartDataPtr<Event::CalDigiCol>(eventSvc(),
                                                                               EventModel::Digi::CalDigiCol);
  if (!calDigiCol) {
    MsgStream msglog(msgSvc(), name());
    msglog << MSG::VERBOSE << "No CalDigi data found"
           << endreq;
    return StatusCode::FAILURE;
  }

  msglog << MSG::DEBUG << "verify CalDigi" << endreq;
  if (m_testCalDigi.verify(*m_calSignalTool,
                           *m_calCalibSvc,
                           m_testCfg,
                           m_twrSet,
                           calDigiCol).isFailure())
    return StatusCode::FAILURE;

  msglog << MSG::DEBUG << "verify CalTrigTool (MC)" << endreq;
  /// use tight tolerance as these answers should be pretty much perfect
  if (m_testCalTrigTool.verify(*m_calSignalTool,
                               m_twrSet,
                               *m_calTrigTool,
                               *m_precalcCalibTool,
                               MAX_RECON_DIFF).isFailure())
    return StatusCode::FAILURE;

  msglog << MSG::DEBUG << "verify CalDiagnosticTool" << endreq;
  if (m_testCalDiagnosticTool.verify(*m_calSignalTool,
                                     *m_calTrigTool,
                                     *m_precalcCalibTool,
                                     *m_calDiagnosticTool,
                                     m_twrSet).isFailure())
    return StatusCode::FAILURE;

  // Retrieve the Event data for this event
  SmartDataPtr<LdfEvent::DiagnosticData> diagTds(eventSvc(), "/Event/Diagnostic");

  msglog << MSG::DEBUG << "verify CalDiagnostic TDS data" << endreq;
  if (m_testCalDiagnosticTool.verifyTDS(diagTds,
                                        *m_calDiagnosticTool,
                                        m_testCfg,
                                        m_twrSet).isFailure())
    return StatusCode::FAILURE;
  
  /// clear out MC so CalTrigToolDigi will work from Digi's instead
  if (clearMC().isFailure())
    return StatusCode::FAILURE;

  msglog << MSG::DEBUG << "verify CalTrigTool (Digi)" << endreq;
  /// use loose tolerance as we are comparing ADC from one diode to
  /// trigger output of another diode
  if (m_testCalTrigTool.verify(*m_calSignalTool,
                               m_twrSet,
                               *m_calTrigToolDigi,
                               *m_precalcCalibTool,
                               MAX_ASYM_DIFF).isFailure())
    return StatusCode::FAILURE;

  // execute CalXtalRecAlg
  sc = m_calXtalRecAlg->execute();
  if (sc.isFailure())
    return sc;

  // retrieve CalXtalRecCol
  // get a pointer to the input TDS data collection
  SmartDataPtr<Event::CalXtalRecCol> calXtalRecCol = SmartDataPtr<Event::CalXtalRecCol>(eventSvc(),
                                                                                        EventModel::CalRecon::CalXtalRecCol);
  if (!calXtalRecCol) {
    MsgStream msglog(msgSvc(), name());
    msglog << MSG::VERBOSE << "No CalXtalRecData found"
           << endreq;
    return StatusCode::FAILURE;
  }

  msglog << MSG::DEBUG << "verify CalXtalRecAlg (MC)" << endreq;
  if (m_testCalXtalRecAlg.verify(*m_calSignalTool,
                                 *m_calCalibSvc,
                                 *m_precalcCalibTool,
                                 m_testCfg,
                                 m_twrSet,
                                 calDigiCol,
                                 calXtalRecCol).isFailure())
    return StatusCode::FAILURE;
                                 
  return StatusCode::SUCCESS;
}

StatusCode test_CalXtalResponse::finalize() {
  return StatusCode::SUCCESS;
}

StatusCode test_CalXtalResponse::applyNewCfg(const TestCfg &testCfg) {
  StatusCode sc;

  // trigMode
  sc = m_calDigiAlg->setProperty("DefaultAllRange", testCfg.trigMode == CalXtalId::ALLRANGE ? "1" : "0");
  if (sc.isFailure())
    return sc;

  // zeroSuppress
  sc = m_calDigiAlg->setProperty("DefaultZeroSuppress", testCfg.zeroSuppress ? "1" : "0");
  if (sc.isFailure())
    return sc;

  // cal diagnostic
  if (m_calDigiAlg->setProperty("CreateDiagnosticData", testCfg.createDiagnosticData ? "1" : "0").isFailure())
    return StatusCode::FAILURE;

  return StatusCode::SUCCESS;
}

/// create mc hit for each hit on each crystal in testCfg
StatusCode test_CalXtalResponse::fillMC(const TestCfg &testCfg, const TwrSet &twrSet) {
  auto_ptr<Event::McIntegratingHitCol> hitCol(new Event::McIntegratingHitCol());

  // xtal loop
  for (TestCfg::XtalList::const_iterator xtalIt(testCfg.testXtals.begin());
       xtalIt != testCfg.testXtals.end();
       xtalIt++) {

    /// don't create mc for crystals not in geometry.
    if (twrSet.find(xtalIt->getTwr()) == twrSet.end())
      continue;

    // hit loop
    for (TestCfg::HitList::const_iterator hitIt(testCfg.xtalHits.begin());
         hitIt != testCfg.xtalHits.end();
         hitIt++) {
      auto_ptr<Event::McIntegratingHit> hit(new Event::McIntegratingHit());

      const float xtalMM = hitIt->xtalPos*m_csiLength;

      fillMcHit(*xtalIt,
                xtalMM,
                hitIt->meV,
                *hit);

      hitCol->push_back(hit.release());
    }
  }

  StatusCode sc = registerMcHitCol(hitCol.release());
  if (sc.isFailure())
    return sc;

  return StatusCode::SUCCESS;
}

StatusCode test_CalXtalResponse::verifyCalCalibSvc(const TestCfg &testCfg) {
  StatusCode  sc;

  /// test 'real' calibrations (currently using partial LAT config to
  /// test proper handling of empty tower bays
  sc = m_testCalCalibSvc.verify(*m_calCalibSvc, m_testCalibSet, testCfg, m_twrSet);
  if (sc.isFailure())
    return sc;

  /// test ideal calibrations which use alternate code path
  sc = m_testCalCalibSvc.verify(*m_calCalibSvcIdeal, m_testCalibSetIdeal, testCfg, m_fullTwrSet);
  if (sc.isFailure())
    return sc;

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
  short nSeg = (short)floor((xtalPos / m_csiLength)             // range 0->1
                            * m_nCsISeg);               // range
                                                        // 0->nCsISeg
            
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
  const float xDiff = xtalPos - segCtr;

  const HepPoint3D mom1(xDiff,0,0);

  //-- Create McIntegratingHit
  hit.setVolumeID(segmId);
  // add single deposit w/ given mev & vector from center 
  // of segment.
  hit.addEnergyItem(meV, NULL, mom1);

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


StatusCode test_CalXtalResponse::clearMC() {
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

  return StatusCode::SUCCESS;
}
