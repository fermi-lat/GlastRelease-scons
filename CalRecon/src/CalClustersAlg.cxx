
#include "CalClustersAlg.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

DECLARE_ALGORITHM_FACTORY(CalClustersAlg) ;

using namespace Event;

CalClustersAlg::CalClustersAlg
 ( const std::string & name, ISvcLocator * pSvcLocator )
 : Algorithm(name,pSvcLocator)
 {   
  // declaration of parameter needed to distinguish 2 calls
  // of CalClustersAlg:
  // 1st - before TkrRecon and 2nd - after TkrRecon
//  declareProperty("callNumber",m_callNumber=0) ;
    
  declareProperty ("clusteringToolName", m_clusteringToolName="CalSingleClustering") ;

  const unsigned int nbDefaultCorrTools = 5 ;
  std::string defaultCorrTools[nbDefaultCorrTools] =
   { "CalLastLayerCorr", "CalProfileCorr", "CalValsCorr",
     "CalTkrLikelihoodCorr", "CalTransvOffsetCorr" } ;
  declareProperty("corrToolNames", m_corrToolNames
   = std::vector<std::string>(defaultCorrTools,defaultCorrTools+nbDefaultCorrTools)) ;

 }



StatusCode CalClustersAlg::initialize()
{
    if (CalReconActor::initialize(serviceLocator()).isFailure()) {
        return StatusCode::FAILURE ;
    }
    
    MsgStream log(msgSvc(), name()) ;
    StatusCode sc = StatusCode::SUCCESS ;
    
    // get callNumber parameter from jobOptions file    
    setProperties();
    
//    log << MSG::INFO << "CalClustersAlg: callNumber = " 
//        << m_callNumber << endreq;
//    
    // set global constants with geometry parameters (in cm)
    // used by Profile() function
    
    sc = toolSvc()->retrieveTool(m_clusteringToolName,  m_clusteringTool);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to create " << m_clusteringToolName << endreq;
        return sc;
    }
        
    const std::vector< std::string > & corrToolNames = m_corrToolNames ;
    std::vector< std::string >::const_iterator toolName ;
    ICalEnergyCorr * tool ;
    for ( toolName = corrToolNames.begin() ;
          toolName != corrToolNames.end() ;
          ++toolName ) {
        sc = toolSvc()->retrieveTool(*toolName,tool) ;
        if (sc.isFailure() ) {
            log<<MSG::ERROR<<" Unable to create "<<*toolName<<endreq ;
            return sc ;
        }
        else {
            m_corrTools.push_back(tool) ;
        }
    }

    return sc;
}

StatusCode CalClustersAlg::execute()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    log<<MSG::DEBUG<<"Begin"<<endreq ;
    getKernel()->reviewEvent() ;
    
    // non fatal errors
    // if there's no CalXtalRec then CalClustersAlg is not happening
  	if (!getKernel()->getXtalRecs())
      return StatusCode::SUCCESS ;
      
    // other errors are fatal
    if (getKernel()->getStatus().isFailure())
      return StatusCode::FAILURE ;
      
    // call the clustering tool
    log<<MSG::DEBUG<<"Clustering"<<endreq ;
    if (m_clusteringTool->findClusters().isFailure())
     {
      log<<MSG::ERROR<<"Clustering failure"<<endreq ;
      return StatusCode::FAILURE ;
     }
    
    // loop over all correction tools
    int itool = 0 ;
    std::vector<ICalEnergyCorr *>::const_iterator tool ;
    for ( tool = m_corrTools.begin() ;
          tool != m_corrTools.end() ;
          ++tool, ++itool )
     {
        log<<MSG::DEBUG<<"Correction "<<itool<<endreq ;
        if (((*tool)->doEnergyCorr()).isFailure())
         { 
          log<<MSG::ERROR<<"Failed correction "<<itool<<endreq ;
          sc = StatusCode::FAILURE ;
         }
     }
        
    // print the reconstruction results for debugging
    getKernel()->getClusters()->writeOut(log<<MSG::DEBUG);     
    
    log<<MSG::DEBUG<<"End"<<endreq ;
    return sc;
}

StatusCode CalClustersAlg::finalize()
 { return StatusCode::SUCCESS ; }




