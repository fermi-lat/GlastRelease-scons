
#include "SingleClusteringTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"

DECLARE_TOOL_FACTORY(SingleClusteringTool) ;

SingleClusteringTool::SingleClusteringTool
 ( const std::string & type, 
   const std::string & name, 
   const IInterface* parent)
 : ClusteringTool(type,name,parent)
 { declareInterface<IClusteringTool>(this) ; }

SingleClusteringTool::~SingleClusteringTool()
 {}

ClusteringTool::xTalDataVec SingleClusteringTool::nextXtalsSet( xTalDataVec & xTalVec )
 {
  xTalDataVec copy = xTalVec ;
  xTalVec.clear() ;
  return copy ;
 }


