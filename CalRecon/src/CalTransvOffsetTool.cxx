
#include "CalTransvOffsetTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_TOOL_FACTORY(CalTransvOffsetTool) ;

CalTransvOffsetTool::CalTransvOffsetTool
 ( const std::string & type,
   const std::string & name,
   const IInterface * parent )
 : EnergyCorr(type,name,parent)
 {
  declareInterface<IEnergyCorr>(this) ;
 }

StatusCode CalTransvOffsetTool::doEnergyCorr( Event::CalCluster * cluster )
 {
  // calculating the transverse offset of average position in the calorimeter
  // with respect to the position predicted from tracker information
  double calTransvOffset = 0.;
  if (getKernel()->getTkrNVertices()>0)
   {
    Vector calOffset = (cluster->getPosition()) - getKernel()->getTkrFrontVertexPosition() ;
    double calLongOffset = getKernel()->getTkrFrontVertexDirection()*calOffset;
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
 

