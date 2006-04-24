// $Header$

// Include files
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "idents/CalXtalId.h"

#include "src/CalFailureModeSvc.h"

#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"
#include "CalUtil/CalArray.h"

#include <vector>

using namespace CalUtil;

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
  //! number of times called
  int m_count; 

  /// pointer to failure mode service
  ICalFailureModeSvc* m_FailSvc;

  StatusCode testCalFailureModeSvc();
  /// test all index & conversion classes in caldefs
  StatusCode testCalDefs();
  StatusCode testCalVec();
  StatusCode testCalArray();
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

  return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode test_CalUtil::execute()
{
  StatusCode  sc = StatusCode::SUCCESS;
  MsgStream   log( msgSvc(), name() );
  log << MSG::INFO << "executing " << ++m_count << " time" << endreq;

  sc = testCalFailureModeSvc();
  if (sc.isFailure()) return sc;
  sc = testCalDefs();
  if (sc.isFailure()) return sc;
  sc = testCalVec();
  if (sc.isFailure()) return sc;
  sc = testCalArray();
  if (sc.isFailure()) return sc;

  
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode test_CalUtil::finalize(){
  StatusCode  sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "finalize after " << m_count << " calls." << endreq;
    
  return sc;
}

StatusCode test_CalUtil::testCalFailureModeSvc() {
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "testCalFailureModeSvc" << endreq;
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

  return StatusCode::SUCCESS;
}

StatusCode test_CalUtil::testCalDefs(){
  MsgStream log(msgSvc(), name());


  log << MSG::INFO << "testCalDefs" << endreq;

  //-- TEST CALDEFS CONVERSION ROUTINES --//
  // try to test all conversions & constructors for
  // consistency.
  for (short tRow=0; tRow < 4 ; tRow++)
    for (short tCol=0; tCol < 4; tCol++) {

      TwrNum twr(tRow, tCol);
      if (tRow != twr.getRow()) return StatusCode::FAILURE;
      if (tCol != twr.getCol()) return StatusCode::FAILURE;
      
      for (DirNum dir; dir.isValid(); dir++)
        for (short dLyr=0; dLyr < 4; dLyr++) {

          LyrNum lyr(dir, dLyr);
          if (dir != lyr.getDir()) return StatusCode::FAILURE;

          for (ColNum col; col.isValid(); col++) {
            CalXtalId xtalId(twr,lyr,col);

            XtalIdx xtalIdx(twr,lyr,col);
            if (xtalId != xtalIdx.getCalXtalId()) return StatusCode::FAILURE;
            if (twr != xtalIdx.getTwr()) return StatusCode::FAILURE;
            if (lyr != xtalIdx.getLyr()) return StatusCode::FAILURE;
            if (col != xtalIdx.getCol()) return StatusCode::FAILURE;

            // alt-constructor(s)
            if (XtalIdx(xtalId) != xtalIdx) return StatusCode::FAILURE;

            for (FaceNum face; face.isValid(); face++) {
              FaceIdx faceIdx(twr,lyr,col,face);
              CalXtalId faceId(twr,lyr,col,(CalXtalId::XtalFace)face);
              if (faceId != faceIdx.getCalXtalId()) return StatusCode::FAILURE;
              if (twr != faceIdx.getTwr()) return StatusCode::FAILURE;
              if (lyr != faceIdx.getLyr()) return StatusCode::FAILURE;
              if (col != faceIdx.getCol()) return StatusCode::FAILURE;
              if (face.val() != faceIdx.getFace().val()) return StatusCode::FAILURE;
              if (xtalIdx != faceIdx.getXtalIdx()) return StatusCode::FAILURE;

              // alt-constructor(s)
              if (FaceIdx(faceId) != faceIdx) return StatusCode::FAILURE;
              if (FaceIdx(xtalIdx, face) != faceIdx) return StatusCode::FAILURE;

              
              for (DiodeNum diode; diode.isValid(); diode++) {
                if (RngNum(diode,THX8).val() != diode.getX8Rng().val()) return StatusCode::FAILURE;
                if (RngNum(diode,THX1) != diode.getX1Rng()) return StatusCode::FAILURE;

                XtalDiode xDiode(face,diode);
                if (diode != xDiode.getDiode()) return StatusCode::FAILURE;
                if (face  != xDiode.getFace()) return StatusCode::FAILURE;

                DiodeIdx diodeIdx(twr,lyr,col,face,diode);
                if (twr != diodeIdx.getTwr()) return StatusCode::FAILURE;
                if (lyr != diodeIdx.getLyr()) return StatusCode::FAILURE;
                if (col != diodeIdx.getCol()) return StatusCode::FAILURE;
                if (face != diodeIdx.getFace()) return StatusCode::FAILURE;
                if (diode != diodeIdx.getDiode()) return StatusCode::FAILURE;
                if (xtalIdx != diodeIdx.getXtalIdx()) return StatusCode::FAILURE;
                if (faceIdx != diodeIdx.getFaceIdx()) return StatusCode::FAILURE;

                //alt-constructor(s)
                if (DiodeIdx(xtalIdx, face, diode) != diodeIdx) return StatusCode::FAILURE;
                if (DiodeIdx(xtalIdx, xDiode) != diodeIdx) return StatusCode::FAILURE;
                if (DiodeIdx(faceIdx, diode) != diodeIdx) return StatusCode::FAILURE;
                
              
                
                for (THXNum thx; thx.isValid(); thx++) {

                  RngNum rng(diode,thx);
                  if (diode != rng.getDiode()) return StatusCode::FAILURE;

                  XtalRng xRng(face, rng);
                  if (face != xRng.getFace()) return StatusCode::FAILURE;
                  if (rng != xRng.getRng()) return StatusCode::FAILURE;
                  if (xDiode != xRng.getXtalDiode()) return StatusCode::FAILURE;

                  RngIdx rngIdx(twr,lyr,col,face,rng);
                  CalXtalId rngId(twr, lyr, col, face.val(), rng.val());
                  if (twr != rngIdx.getTwr()) return StatusCode::FAILURE;
                  if (lyr != rngIdx.getLyr()) return StatusCode::FAILURE;
                  if (col != rngIdx.getCol()) return StatusCode::FAILURE;
                  if (face != rngIdx.getFace()) return StatusCode::FAILURE;
                  if (rng != rngIdx.getRng()) return StatusCode::FAILURE;
                  if (xtalIdx != rngIdx.getXtalIdx()) return StatusCode::FAILURE;
                  if (faceIdx != rngIdx.getFaceIdx()) return StatusCode::FAILURE;
                  
                  //alt-constructor(s)
                  if (RngIdx(xtalIdx, face, rng) != rngIdx) return StatusCode::FAILURE;
                  if (RngIdx(xtalIdx, xRng) != rngIdx) return StatusCode::FAILURE;
                  if (RngIdx(faceIdx, rng) != rngIdx) return StatusCode::FAILURE;


                } // thx/rng
              } // diode
            } // face
          } // col
        } // dLyr
    } // twr

  return StatusCode::SUCCESS;
}

StatusCode test_CalUtil::testCalVec(){
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "testCalVec" << endreq;

  // creats vector of size 8, which is
  // enought to test functionality
  CalVec<XtalRng, unsigned> testVec;

  testVec.resize(XtalRng::N_VALS);
  
  if (testVec.size() != (unsigned)XtalRng::N_VALS)
    return StatusCode::FAILURE;

  for (XtalRng idx; idx.isValid(); idx++) {
    testVec[idx] = idx.val();
  }

  for (XtalRng idx; idx.isValid(); idx++) {
    if (testVec[idx] != idx.val())
      return StatusCode::FAILURE;
    
    if (testVec.find(idx.val()) - testVec.begin() 
        != (int)idx.val())
      return StatusCode::FAILURE;
  }

  testVec.fill(3);
  for (XtalRng idx; idx.isValid(); idx++) {
    if (testVec[idx] != 3)
      return StatusCode::FAILURE;
  }

  testVec.clear();
  if (testVec.size() != (unsigned)0) 
    return StatusCode::FAILURE;

  return StatusCode::SUCCESS;
}

StatusCode test_CalUtil::testCalArray(){
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "testCalArray" << endreq;

  // creats array of size 8, which is
  // enough to test functionality
  CalArray<XtalRng, unsigned> testArr;

  if (testArr.size() != (unsigned)XtalRng::N_VALS)
    return StatusCode::FAILURE;

  for (XtalRng idx; idx.isValid(); idx++) {
    testArr[idx] = idx.val();
  }

  for (XtalRng idx; idx.isValid(); idx++) {
    if (testArr[idx] != idx.val())
      return StatusCode::FAILURE;
    
    if ((testArr.find(idx.val()) - testArr.begin()) 
        != (int)idx.val())
      return StatusCode::FAILURE;
  }

  testArr.fill(3);
  for (XtalRng idx; idx.isValid(); idx++) {
    if (testArr[idx] != 3)
      return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}
