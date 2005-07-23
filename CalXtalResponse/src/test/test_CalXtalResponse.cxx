// $Header$

// Include files
// Gaudi system includes

// LOCAL
#include "CalXtalResponse/ICalCalibSvc.h"
#include "CalXtalResponse/IXtalRecTool.h"
#include "CalXtalResponse/IXtalDigiTool.h"

// GLAST
#include "idents/CalXtalId.h"
#include "CalUtil/CalDefs.h"

// EXTLIB
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

// STD

using namespace CalDefs;

// Define the class here instead of in a header file: not needed anywhere but here!
//------------------------------------------------------------------------------
/** 
    A simple algorithm.

*/
class test_CalXtalResponse : public Algorithm {
public:
  test_CalXtalResponse(const string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

private: 
  StatusCode test_calCalibSvc();
  StatusCode test_xtalRecTool();
  StatusCode test_xtalDigiTool();

  //! number of times called
  int m_count; 

  IXtalRecTool *m_xtalRecTool; ///< pointer to calenergyTool
  IXtalDigiTool *m_xtalDigiTool; ///< pointer to XtalDigiTool

};
//------------------------------------------------------------------------

// necessary to define a Factory for this algorithm
// expect that the xxx_load.cxx file contains a call     
//     DLL_DECL_ALGORITHM( test_CalXtalResponse );

static const AlgFactory<test_CalXtalResponse>  Factory;
const IAlgFactory& test_CalXtalResponseFactory = Factory;

//------------------------------------------------------------------------
//! ctor
test_CalXtalResponse::test_CalXtalResponse(const string& name, ISvcLocator* pSvcLocator)
  :Algorithm(name, pSvcLocator)
  ,m_count(0)
{
}

//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode test_CalXtalResponse::initialize(){
  StatusCode  sc;
  MsgStream msglog(msgSvc(), name()); 
  msglog << MSG::INFO << "initialize" << endreq;

  // get cal energy Tool
  sc = toolSvc()->retrieveTool("XtalRecTool", m_xtalRecTool);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "XtalRecTool" << endreq;
    return sc;
  }

  // get cal adc Tool
  sc = toolSvc()->retrieveTool("XtalDigiTool", m_xtalDigiTool);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << "XtalDigiTool" << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------
//! process an event
StatusCode test_CalXtalResponse::execute()
{
  MsgStream msglog(msgSvc(), name());
  msglog << MSG::INFO << "executing " << ++m_count << " time" << endreq;
  StatusCode sc;

  // Run each individual test
  if ((sc = test_calCalibSvc().isFailure())) return sc;
  if ((sc = test_xtalRecTool().isFailure())) return sc;
  if ((sc = test_xtalDigiTool().isFailure())) return sc;

  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode test_CalXtalResponse::finalize(){
  MsgStream msglog(msgSvc(), name());
  msglog << MSG::INFO << "finalize after " << m_count << " calls." << endreq;

  return StatusCode::SUCCESS;
}

/// test one of each calib_type from CalCalibSvc
StatusCode test_CalXtalResponse::test_calCalibSvc() {
  MsgStream         msglog( msgSvc(), name() );
  StatusCode sc;

  // Setup CalCalibSvc
  ICalCalibSvc *pCalCalibSvc;
  sc = service("CalCalibSvc", pCalCalibSvc, false);
  if (!pCalCalibSvc) {
    msglog << MSG::ERROR << "can't get CalCalibSvc" << endreq;
    return StatusCode::FAILURE;
  }

  // pick one xtal/rng combo
  int lyr = 2;
  int col = 3;
  FaceNum face = POS_FACE;
  RngNum rng = LEX1;
  CalXtalId xtalId(0,lyr,col,face,rng);

  msglog << MSG::INFO << "Testing CalCalibSvc lyr=" << lyr 
         << " col="  << col 
         << " face=" << face
         << " rng="  << rng 
         << endreq;

  // IntNonlin
  const vector<float> *vals;
  const vector<unsigned> *dacs;
  float error;
  if ((sc = pCalCalibSvc->getIntNonlin(xtalId, 
                                       vals, dacs, error)).isFailure()) {
    msglog << MSG::ERROR << "Error retrieving intNonlin from CalCalibSvc" << endreq;
    return sc;
  }    
  // if we gotem, then let's loop through
  for (unsigned i = 0; i < vals->size(); i++)
    msglog << MSG::INFO << "Calib: CalXtalId=" << xtalId
           << " Dac=" << (*dacs)[i]
           << " ADC=" << (*vals)[i]
           << endreq;

  // PED
  float ped, cos, sig;
  if ((sc = pCalCalibSvc->getPed(xtalId, 
                                 ped, sig, cos)).isFailure()) {
    msglog << MSG::ERROR << "Error retrieving ped from CalCalibSvc" << endreq;
    return sc;
  }
  msglog << MSG::INFO << "Calib: CalXtalId=" << xtalId
         << " ped="  << ped
         << " sig=" << sig
         << endreq;

  // ASYM
  CalXtalId xtalId2(0, lyr, col);
  const vector<CalibData::ValSig> *lrgVec, *smVec, *nspbVec, *psnbVec;
  const vector<float> *xVals;
  if ((sc = pCalCalibSvc->getAsym(xtalId2, lrgVec,smVec,nspbVec,psnbVec,xVals)) 
      .isFailure()) {
    msglog << MSG::ERROR << "Error retrieving Asym from CalCalibSvc" << endreq;
    return sc;
  }
  msglog << MSG::DEBUG << "Asymmetry successfully getd" << endreq;

  // MeVPerDac
  CalibData::ValSig mpdLrg, mpdSm;
  if ((sc = pCalCalibSvc->getMeVPerDac(xtalId2, mpdLrg, mpdSm)).isFailure()) {
    msglog << MSG::ERROR << "Error retrieving MPD from CalCalibSvc" << endreq;
    return sc;
  }
  msglog << MSG::DEBUG << "MeV per Dac successfully getd" << endreq;

  return StatusCode::SUCCESS;
}

/// do single digi->energy calculation
StatusCode test_CalXtalResponse::test_xtalRecTool() {
#if 0
  MsgStream         msglog( msgSvc(), name() );
  StatusCode sc;

  // pick one xtal/rng combo
  CalXtalId::XtalFace face = CalXtalId::POS;
  CalXtalId::AdcRange rng = CalXtalId::LEX1;
  CalXtalId xtalId(0,2,3,face,rng);

  bool belowThreshP, belowThreshN, saturatedP, saturatedN;
  double pos, energy;
  if ((sc = m_xtalRecTool->calculate(xtalId, rng, rng, 
                                     1000, 1000, energy, pos,
                                     belowThreshP,
                                     belowThreshN,
                                     saturatedP,
                                     saturatedN))
      != sc) {
    msglog << MSG::ERROR << "Error calling XtalRecTool::calculate()" << endreq;
    return sc;
  }
  msglog << MSG::INFO << "2 face energy calculated:" << energy << endreq;
  bool belowThresh, saturated;
  if ((sc = m_xtalRecTool->calcSingleFaceEne(xtalId, 1000, 0, energy, 
                                     belowThresh, saturated))
      != sc) {
    msglog << MSG::ERROR << "Error calling XtalRecTool::calculate()" << endreq;
    return sc;
  }
  msglog << MSG::INFO << "1 face energy calculated:" << energy << endreq;
#endif

  return StatusCode::SUCCESS;
}

/// do single mc->digi calculation, 
/// I'm too lazy to build my own hits, so i'm just sending empty list for now.
StatusCode test_CalXtalResponse::test_xtalDigiTool() {
#if 0
  MsgStream msglog(msgSvc(), name());
  StatusCode sc;

  CalXtalId xtalId(0,2,3);

  CalXtalId::AdcRange rngP, rngN;
  bool lacP, lacN;
  vector<int> adcP(4,0), adcN(4,0);

  // call XtalDigiTool
  vector<const Event::McIntegratingHit*> hitvec(0);
  
  

  if ((sc = m_xtalDigiTool->calculate(xtalId,
                                     hitvec,
                                     lacP, lacN,
                                     rngP, rngN,
                                     adcP, adcN
                                     )).isFailure()) {
    msglog << MSG::ERROR << "error calling XtalDigiTool" << endreq;
    return sc;
  }

  msglog << MSG::INFO << "XtalDigiTool " << xtalId
         << " lacP " << lacP
         << " rngP " << rngP
         << " adcP " << adcP[0]
         << ' '      << adcP[1]
         << ' '      << adcP[2]
         << ' '      << adcP[3]
         << " lacN " << lacN
         << " rngN " << rngN
         << " adcN " << adcN[0]
         << ' '      << adcN[1]
         << ' '      << adcN[2]
         << ' '      << adcN[3]
         << endreq;
#endif
  return StatusCode::SUCCESS;
}
