
#include "CalTransvOffsetCorr.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_TOOL_FACTORY(CalTransvOffsetCorr) ;

CalTransvOffsetCorr::CalTransvOffsetCorr
 ( const std::string & type,
   const std::string & name,
   const IInterface * parent )
 : CalEnergyCorr(type,name,parent)
 {
  declareInterface<ICalEnergyCorr>(this) ;
 }

StatusCode CalTransvOffsetCorr::doEnergyCorr( Event::CalCluster * cluster )
 {
  // calculating the transverse offset of average position in the calorimeter
  // with respect to the position predicted from tracker information
  double calTransvOffset = 0.;
  if (getKernel()->getTkrNVertices()>0)
   {
    Vector calOffset = (cluster->getPosition()) - getKernel()->getTkrFrontVertex()->getPosition() ;
    double calLongOffset = getKernel()->getTkrFrontVertex()->getDirection()*calOffset;
    calTransvOffset = sqrt(calOffset.mag2() - calLongOffset*calLongOffset);
   }
    
  // store the calculated quantities back in this CalCluster object. Note
  // temporary ugly kluge of overwriting all but two existing elements!
  cluster->initialize(
    cluster->getEnergyLeak(),
    cluster->getEneLayer(),
    cluster->getPosLayer(),
    cluster->getRmsLayer(),
    cluster->getRmsLong(),
    cluster->getRmsTrans(),
    cluster->getDirection(),
    calTransvOffset) ;
    
  return StatusCode::SUCCESS ;
 }
 

