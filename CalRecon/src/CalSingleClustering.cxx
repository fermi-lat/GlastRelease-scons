
#include "CalSingleClustering.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"

DECLARE_TOOL_FACTORY(CalSingleClustering) ;

CalSingleClustering::CalSingleClustering
 ( const std::string & type, 
   const std::string & name, 
   const IInterface* parent)
 : CalClustering(type,name,parent)
 { declareInterface<ICalClustering>(this) ; }

CalSingleClustering::~CalSingleClustering()
 {}

void CalSingleClustering::makeSets( const XtalDataVec & xtals, XtalDataVecVec & clusters )
 { clusters.push_back(new XtalDataVec(xtals)) ; }


