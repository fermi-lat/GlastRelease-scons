// $Header$

// Include files
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "idents/CalXtalId.h"

#include "CalUtil/CalFailureModeSvc.h"
#include "CalUtil/ICalCalibSvc.h"
#include "CalUtil/ICalEnergyTool.h"
#include "CalUtil/ICalPosTool.h"
#include "CalUtil/ICalAdcTool.h"

// Define the class here instead of in a header file: not needed anywhere but here!
//------------------------------------------------------------------------------
/** 
    A simple algorithm.

  
*/
class test_CalUtil : public Algorithm {
public:
  test_CalUtil(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();
    
private: 
  StatusCode test_calCalibSvc();
  StatusCode test_testEnergyTool();
  StatusCode test_testPosTool();
  StatusCode test_testAdcTool();

  //! number of times called
  int m_count; 

  /// pointer to failure mode service
  ICalFailureModeSvc* m_FailSvc;

  ICalEnergyTool *m_pTestEnergyTool; ///< pointer to testenergytool
  ICalPosTool *m_pTestPosTool; ///< pointer to testpostool
  ICalAdcTool *m_pTestAdcTool; ///< pointer to testadctool

};
//------------------------------------------------------------------------

// necessary to define a Factory for this algorithm
// expect that the xxx_load.cxx file contains a call     
//     DLL_DECL_ALGORITHM( test_CalUtil );

static const AlgFactory<test_CalUtil>  Factory;
const IAlgFactory& test_CalUtilFactory = Factory;

//------------------------------------------------------------------------
//! ctor
test_CalUtil::test_CalUtil(const std::string& name, ISvcLocator* pSvcLocator)
  :Algorithm(name, pSvcLocator)
  ,m_count(0)
{
}

//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode test_CalUtil::initialize(){
  StatusCode  sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "initialize" << endreq;
    
  sc = service("CalFailureModeSvc", m_FailSvc);
  if (sc.isFailure() ) {
    log << MSG::ERROR << "  Unable to find CalFailureMode service" << endreq;
    return sc;
  }

  // get cal energy tool
  sc = toolSvc()->retrieveTool("TestEnergyTool", m_pTestEnergyTool);
  if (sc.isFailure() ) {
    log << MSG::ERROR << "  Unable to create " << "TestEnergyTool" << endreq;
    return sc;
  }

  // get cal pos tool
  sc = toolSvc()->retrieveTool("TestPosTool", m_pTestPosTool);
  if (sc.isFailure() ) {
    log << MSG::ERROR << "  Unable to create " << "TestPosTool" << endreq;
    return sc;
  }

  // get cal adc tool
  sc = toolSvc()->retrieveTool("TestAdcTool", m_pTestAdcTool);
  if (sc.isFailure() ) {
    log << MSG::ERROR << "  Unable to create " << "TestAdcTool" << endreq;
    return sc;
  }

  return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode test_CalUtil::execute()
{
  StatusCode  sc = StatusCode::SUCCESS;
  MsgStream   log( msgSvc(), name() );
  log << MSG::INFO << "executing " << ++m_count << " time" << endreq;

  idents::CalXtalId id1(10,2,2);  // tower 10, layer 3(y)
  idents::CalXtalId id2(11,1,2);  // tower 11, layer 1(y)
  idents::CalXtalId id3(3,5,3);   // tower 3, layer 5(y)
  idents::CalXtalId id4(4,6,3);   // tower 4, layer 6(x)
  idents::CalXtalId id5(4,1,3);   // tower 4, layer 1(y)

  if (m_FailSvc == 0) return StatusCode::FAILURE;
  if (m_FailSvc->matchChannel(id1,idents::CalXtalId::NEG)) {
    log << MSG::INFO << "removed channel (10,2,2) NEG" << endreq;
  } else {
    log << MSG::ERROR << "failed to remove channel (10,2,2) NEG" << endreq;
    return StatusCode::FAILURE;
  }
  if (m_FailSvc->matchChannel(id2,idents::CalXtalId::NEG)) {
    log << MSG::INFO << "removed channel (11,1,2) NEG" << endreq;
  } else {
    log << MSG::ERROR << "failed to remove channel (11,1,2) NEG" << endreq;
    return StatusCode::FAILURE;
  }
  if (!m_FailSvc->matchChannel(id3,idents::CalXtalId::NEG)) {
    log << MSG::INFO << "left channel (3,5,3) NEG" << endreq;
  } else {
    log << MSG::ERROR << "erroneously removed channel (3,5,3) NEG" << endreq;
    return StatusCode::FAILURE;
  }
  if (m_FailSvc->matchChannel(id4,idents::CalXtalId::NEG)) {
    log << MSG::INFO << "removed channel (4,6,3) NEG" << endreq;
  } else {
    log << MSG::ERROR << "erroneously left channel (4,6,3) NEG" << endreq;
    return StatusCode::FAILURE;
  }
  if (!m_FailSvc->matchChannel(id5,idents::CalXtalId::NEG)) {
    log << MSG::INFO << "left channel (4,1,3) NEG" << endreq;
  } else {
    log << MSG::ERROR << "erroneously removed channel (4,1,3) NEG" << endreq;
    return StatusCode::FAILURE;
  }

  if (test_calCalibSvc() != StatusCode::SUCCESS) return StatusCode::FAILURE;
  if (test_testEnergyTool() != StatusCode::SUCCESS) return StatusCode::FAILURE;
  if (test_testPosTool() != StatusCode::SUCCESS) return StatusCode::FAILURE;
  if (test_testAdcTool() != StatusCode::SUCCESS) return StatusCode::FAILURE;

  return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode test_CalUtil::finalize(){
  StatusCode  sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "finalize after " << m_count << " calls." << endreq;
    
  return sc;
}

/// test one of each calib_type from CalCalibSvc
StatusCode test_CalUtil::test_calCalibSvc() {
  MsgStream         log( msgSvc(), name() );
  StatusCode sc = StatusCode::SUCCESS;
  
  // Setup CalCalibSvc
  ICalCalibSvc *pCalCalibSvc;
  sc = service("CalCalibSvc", pCalCalibSvc, false);
  if (!pCalCalibSvc) {
    log << MSG::ERROR << "Unable to retrieve CalCalibSvc" << endreq;
    return StatusCode::FAILURE;
  }

  // pick one xtal/range combo
  idents::CalXtalId::XtalFace face = idents::CalXtalId::POS;
  idents::CalXtalId::AdcRange range = idents::CalXtalId::LEX1;
  idents::CalXtalId xtalId(0,2,3);

  const std::vector<float> *vals;
  const std::vector<unsigned> *dacs;
  float error;
  if (pCalCalibSvc->getIntNonlin(xtalId, 
                                 face,
                                 range, 
                                 vals, dacs, error) != StatusCode::SUCCESS) {
    log << MSG::ERROR << "Error retrieving intNonlin from CalCalibSvc" << endreq;
    return StatusCode::FAILURE;
  }
                  
  // if we gotem, then let's loop through
  for (unsigned i = 0; i < vals->size(); i++) {
    log << MSG::INFO << "CalCalib: xtalid=" << xtalId 
        << " face="  << face
        << " range=" << range
        << " Dac="   << (*dacs)[i]
        << " ADC="   << (*vals)[i]
        << endreq;

  }

  // GAIN
  float gain, sig;
  if (pCalCalibSvc->getGain(xtalId, 
                            face,
                            range, 
                            gain, sig) != StatusCode::SUCCESS) {
    log << MSG::ERROR << "Error retrieving gain from CalCalibSvc" << endreq;
    return StatusCode::FAILURE;
  } else {
    log << MSG::INFO << "CalCalib: xtalid=" << xtalId 
        << " face="  << face
        << " range=" << range
        << " gain="  << gain
        << " error=" << error
        << endreq;
  }

  // LightAsym
  const std::vector<float> *lightAsym;
  if (pCalCalibSvc->getLightAsym(xtalId, 
                                 face,
                                 range, 
                                 lightAsym, sig) != StatusCode::SUCCESS) {
    log << MSG::ERROR << "Error retrieving LightAsym from CalCalibSvc" << endreq;
    return StatusCode::FAILURE;
  } else {
    log << MSG::INFO << "Found " << lightAsym->size() << " light asymmetry values." << endreq;
    for (unsigned i = 0; i < lightAsym->size(); i++) {
      log << MSG::INFO << "CalCalib: xtalid=" << xtalId 
          << " face="  << face
          << " range=" << range
          << " LightAsym[" << i << "]="  << (*lightAsym)[i]
          << " error=" << error
          << endreq;
    }
  }

  // MUSLOPE
  float muSlope;
  if (pCalCalibSvc->getMuSlope(xtalId, 
                               face,
                               range, 
                               muSlope, sig) != StatusCode::SUCCESS) {
    log << MSG::ERROR << "Error retrieving muSlope from CalCalibSvc" << endreq;
    return StatusCode::FAILURE;
  } else {
    log << MSG::INFO << "CalCalib: xtalid=" << xtalId 
        << " face="  << face
        << " range=" << range
        << " muSlope="  << muSlope
        << " error=" << error
        << endreq;
  }

  // PED
  float ped, cos;
  if (pCalCalibSvc->getPed(xtalId, 
                           face,
                           range, 
                           ped, sig, cos) != StatusCode::SUCCESS) {
    log << MSG::ERROR << "Error retrieving ped from CalCalibSvc" << endreq;
    return StatusCode::FAILURE;
  } else {
    log << MSG::INFO << "CalCalib: xtalid=" << xtalId 
        << " face="  << face
        << " range=" << range
        << " ped="  << ped
        << " sig=" << sig
        << endreq;
  }

  // GET ROOT TSPLINE for intNonlin
  const TSpline3 *intNonlinSpline;
  if (pCalCalibSvc->getIntNonlin(xtalId, 
                                 face,
                                 range, 
                                 intNonlinSpline) != StatusCode::SUCCESS)
    log << MSG::ERROR << "Error retrieving ped from CalCalibSvc" << endreq;
  //else intNonlinSpline->Dump();

  return StatusCode::SUCCESS;
}

/// do single digi->energy calculation
StatusCode test_CalUtil::test_testEnergyTool() {
  MsgStream         log( msgSvc(), name() );

  // pick one xtal/range combo
  idents::CalXtalId::XtalFace face = idents::CalXtalId::POS;
  idents::CalXtalId::AdcRange range = idents::CalXtalId::LEX1;
  idents::CalXtalId xtalId(0,2,3);
  float energy;

  if (m_pTestEnergyTool->calculate(xtalId, range, range, 200, 200, .5, energy)
      != StatusCode::SUCCESS) {
    log << MSG::ERROR << "Error calling TestEnergyTool::calculate()" << endreq;
    return StatusCode::FAILURE;
  }
  log << MSG::INFO << "2 face energy calculated:" << energy << endreq;

  if (m_pTestEnergyTool->calculate(xtalId, face, range, 200, .5, energy)
      != StatusCode::SUCCESS) {
    log << MSG::ERROR << "Error calling TestEnergyTool::calculate()" << endreq;
    return StatusCode::FAILURE;
  }
  log << MSG::INFO << "1 face energy calculated:" << energy << endreq;

  return StatusCode::SUCCESS;
}

/// do single digi->pos calculation
StatusCode test_CalUtil::test_testPosTool() {
  MsgStream log(msgSvc(), name());

  // pick one xtal/range combo
  idents::CalXtalId::AdcRange range = idents::CalXtalId::LEX1;
  idents::CalXtalId xtalId(0,2,3);
  float pos;

  if (m_pTestPosTool->calculate(xtalId, range, range, 200, 250, pos)
      != StatusCode::SUCCESS) {
    log << MSG::ERROR << "Error calling TestPosTool::calculate()" << endreq;
    return StatusCode::FAILURE;
  }
  log << MSG::INFO << "position calculated:" << pos << endreq;
  
  return StatusCode::SUCCESS;
}

/// do single mc->digi calculation, 
/// I'm too lazy to build my own hits, so i'm just sending empty list for now.
StatusCode test_CalUtil::test_testAdcTool() {
  MsgStream log(msgSvc(), name());

  idents::CalXtalId xtalId(0,2,3);

  idents::CalXtalId::AdcRange rangeP, rangeN;
  bool lacP, lacN;
  std::vector<int> adcP(4,0), adcN(4,0);
  
  
  // call CalAdcTool
  std::vector<const Event::McIntegratingHit*> hitvec(0);
  
  if (m_pTestAdcTool->calculate(xtalId,
                                hitvec,
                                lacP, lacN,
                                rangeP, rangeN,
                                adcP, adcN
                                ) != StatusCode::SUCCESS) {
    log << MSG::ERROR << "error calling testAdcTool" << endreq;
    return StatusCode::FAILURE;
  }
  
  log << MSG::INFO << "CalAdcTool " << xtalId
      << " lacP " << lacP
      << " rangeP " << rangeP
      << " adcP " << adcP[0]
      << " "     << adcP[1]
      << " "     << adcP[2]
      << " "     << adcP[3]
      << " lacN " << lacN
      << " rangeN " << rangeN
      << " adcN " << adcN[0]
      << " "     << adcN[1]
      << " "     << adcN[2]
      << " "     << adcN[3]
      << endreq;

  return StatusCode::SUCCESS;
}
