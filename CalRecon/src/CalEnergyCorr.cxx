
#include "CalEnergyCorr.h"

StatusCode CalEnergyCorr::initialize()
 {
  MsgStream log(msgSvc(),name()) ;
  if (service("CalReconSvc",m_calReconSvc,true).isFailure())
   {
    log<<MSG::ERROR<<"Could not find CalReconSvc"<<endreq ;
    return StatusCode::FAILURE ;
   }
  return StatusCode::SUCCESS ;
 }

StatusCode CalEnergyCorr::doEnergyCorr()
 {
  StatusCode sc = StatusCode::SUCCESS ;
  int icluster ;
  Event::CalClusterCol::const_iterator it ;
  for ( it = m_calReconSvc->getClusters()->begin(), icluster=0 ;
        it != m_calReconSvc->getClusters()->end() ;
        ++it, ++icluster )
   {
    if (doEnergyCorr(*it).isFailure())
     {
      MsgStream log(msgSvc(),name()) ;
      log<<MSG::ERROR<<"Failed correction on cluster "<<icluster<<endreq ;
      sc = StatusCode::FAILURE ;
     }
   }
  return sc ;
 }

