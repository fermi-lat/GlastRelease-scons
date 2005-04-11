
#include "EnergyCorr.h"

StatusCode EnergyCorr::initialize()
 {
  if (CalReconActor::initialize(serviceLocator()).isFailure())
   { return StatusCode::FAILURE ; }
  return StatusCode::SUCCESS ;
 }

