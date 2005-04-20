
#include "CalEnergyCorr.h"

StatusCode CalEnergyCorr::initialize()
 {
  if (CalReconActor::initialize(serviceLocator()).isFailure())
   { return StatusCode::FAILURE ; }
  return StatusCode::SUCCESS ;
 }

StatusCode CalEnergyCorr::doEnergyCorr()
 {
  StatusCode sc = StatusCode::SUCCESS ;
  int icluster ;
  Event::CalClusterCol::const_iterator it ;
  for ( it = getKernel()->getClusters()->begin(), icluster=0 ;
        it != getKernel()->getClusters()->end() ;
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

