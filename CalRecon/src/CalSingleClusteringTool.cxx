
#include "CalSingleClusteringTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"

DECLARE_TOOL_FACTORY(CalSingleClusteringTool) ;

CalSingleClusteringTool::CalSingleClusteringTool
 ( const std::string & type, 
   const std::string & name, 
   const IInterface* parent)
 : CalClusteringTool(type,name,parent)
 { declareInterface<CalIClusteringTool>(this) ; }

CalSingleClusteringTool::~CalSingleClusteringTool()
 {}

void CalSingleClusteringTool::makeSets( const XtalDataVec & xtals, XtalDataVecVec & clusters )
 { clusters.push_back(new XtalDataVec(xtals)) ; }


