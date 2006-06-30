//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Anc/AncCalibTaggerGain.h"
#include "CalibData/Anc/AncCalibTaggerPed.h"
#include "CalibData/Anc/AncCalibQdcPed.h"

/**
   @file UseAnc.cxx                         
   Simple algorithm to test functioning of "the other" TDS with
   ancillary detector calib. data
*/



  /** 
   @class UseAnc

   Algorithm exemplifying retrieval and use of ancillary calibrations
*/
class UseAnc : public Algorithm {


public:
  UseAnc(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Helper functions called by execute
  void processNew(CalibData::AncCalibTaggerPed* pNewTagger,
                  const std::string& pathTagger,
                  CalibData::AncCalibQdcPed* pNewQdc, 
                  const std::string& pathQdc);

  IDataProviderSvc* m_pCalibDataSvc;
  int               m_serTaggerPed;
  int               m_serTaggerGain;   // unused for now
  int               m_serQdcPed;
};


/// Instantiation of a static factory to create instances of this algorithm
static const AlgFactory<UseAnc> Factory;
const IAlgFactory& UseAncFactory = Factory;


UseAnc::UseAnc(const std::string&  name, 
                 ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
    m_serTaggerPed(-1), m_serTaggerGain(-1),m_serQdcPed(-1)
{
  // Declare properties here.

}


StatusCode UseAnc::initialize() {
  StatusCode sc;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initialize()" << endreq;

  // So far don't have any properties, but in case we do some day..
  setProperties();


  sc = service("CalibDataSvc", m_pCalibDataSvc, true);

  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
	<< "Could not get IDataProviderSvc interface of CalibDataSvc" 
	<< endreq;
    return sc;
  }

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  return StatusCode::SUCCESS;

}


StatusCode UseAnc::execute( ) {

  MsgStream log(msgSvc(), name());

  std::string fullPathTaggerPed = "/Calib/ANC_TaggerPed/vanilla";
  std::string fullPathQdcPed = "/Calib/ANC_QdcPed/vanilla";
  DataObject *pObjectTagger;
  DataObject *pObjectQdc;
  

  m_pCalibDataSvc->retrieveObject(fullPathTaggerPed, pObjectTagger);
  m_pCalibDataSvc->retrieveObject(fullPathQdcPed, pObjectQdc);

  CalibData::AncCalibTaggerPed* pTaggerPeds = 0;
  CalibData::AncCalibQdcPed* pQdcPeds = 0;
  pTaggerPeds = dynamic_cast<CalibData::AncCalibTaggerPed *> (pObjectTagger);
  if (!pTaggerPeds) {
    log << MSG::ERROR << "Dynamic cast to AncCalibTaggerPed failed" << endreq;
    return StatusCode::FAILURE;
  }

  pQdcPeds = dynamic_cast<CalibData::AncCalibQdcPed *> (pObjectQdc);
  if (!pQdcPeds) {
    log << MSG::ERROR << "Dynamic cast to AncCalibQdcPed failed" << endreq;
    return StatusCode::FAILURE;
  }

  int newSerNoTagger = pTaggerPeds->getSerNo();
  int newSerNoQdc = pQdcPeds->getSerNo();
  if ((newSerNoTagger != m_serTaggerPed) ||
      (newSerNoQdc != m_serQdcPed) )  {
    log << MSG::INFO << "Processing new gains, peds after retrieveObject" 
        << endreq;
    m_serTaggerPed = newSerNoTagger;
    m_serQdcPed = newSerNoQdc;
    processNew(pTaggerPeds, fullPathTaggerPed, pQdcPeds, fullPathQdcPed);
  }
  m_pCalibDataSvc->updateObject(pObjectTagger);
  m_pCalibDataSvc->updateObject(pObjectQdc);


  pTaggerPeds = 0;
  pQdcPeds = 0;
  try {
    pTaggerPeds = 
      dynamic_cast<CalibData::AncCalibTaggerPed *> (pObjectTagger);
    pQdcPeds = dynamic_cast<CalibData::AncCalibQdcPed *> (pObjectQdc);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to AncCalibTaggerPed or AncCalibQdcPed after " 
        << "update failed" 
        << endreq;
    return StatusCode::FAILURE;
  }
  newSerNoTagger = pTaggerPeds->getSerNo();
  newSerNoQdc = pQdcPeds->getSerNo();
  if ((newSerNoTagger != m_serTaggerPed) ||
      (newSerNoQdc != m_serQdcPed) )  {
    log << MSG::INFO << "Processing new tagger & qdc peds after update" 
        << endreq;
    m_serTaggerPed = newSerNoTagger;
    m_serQdcPed = newSerNoQdc;
    processNew(pTaggerPeds, fullPathTaggerPed, pQdcPeds, fullPathQdcPed);
  }

  return StatusCode::SUCCESS;
}


void UseAnc::processNew(CalibData::AncCalibTaggerPed* pTagger, 
                        const std::string& pathTagger,
                        CalibData::AncCalibQdcPed* pQdc,
                        const std::string& pathQdc) {
  using CalibData::AncTaggerPed;
  using CalibData::AncQdcPed;
  bool  done = false;

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Retrieved with path " << pathTagger << endreq
      << "Serial #" <<  pTagger->getSerNo() << endreq; 
  log << MSG::INFO << "Vstart: " <<  (pTagger->validSince()).hours()
      << "  Vend: " << (pTagger->validTill()).hours() << endreq;
  
  log << MSG::INFO << "And retrieved with path " << pathQdc << endreq
      << "Serial #" <<  pQdc->getSerNo() << endreq; 
  log << MSG::INFO << "Vstart: " <<  (pQdc->validSince()).hours()
      << "  Vend: " << (pQdc->validTill()).hours() << endreq;
  
  if (!done) {
    done = true;
    int iMod = 0;
    int iLay = 0;
    int iChan = 1;

    CalibData::RangeBase* pTaggerPed = pTagger->getChan(iMod, iLay, iChan);
    CalibData::RangeBase* pQdcPed = pQdc->getQdcChan(iMod, iChan);
    
    AncTaggerPed* pT = dynamic_cast<AncTaggerPed * >(pTaggerPed);
    AncQdcPed* pQ = dynamic_cast<AncQdcPed * >(pQdcPed);

    log << MSG::INFO << "For tagger module = " << iMod << " layer = " << iLay
        << " channel = " << iChan << endreq;
    
    log << MSG::INFO << "ped = " << pT->getVal() << endreq;
    log << MSG::INFO << "rNoise = " << pT->getRNoise() << endreq;
    log << MSG::INFO << "sNoise = " << pT->getSNoise() << endreq;
    log << MSG::INFO << "isBad = " << pT->getIsBad() << endreq << endreq;


    log << MSG::INFO << "For qdc module = " << iMod 
        << " channel = " << iChan << endreq;
    log << MSG::INFO << " ped = " << pQ->getVal() << endreq;
    log << MSG::INFO << " rms = " << pQ->getRms() << endreq;
    log << MSG::INFO << " isBad = " << pQ->getIsBad() << endreq;
    log << MSG::INFO << " device = " << pQ->getDevice() << endreq;
  }
}

StatusCode UseAnc::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "          Finalize UseAnc "
      << endreq;
  
  return StatusCode::SUCCESS;
}

