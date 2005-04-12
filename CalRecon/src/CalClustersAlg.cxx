
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
    
  declareProperty ("clusteringToolName", m_clusteringToolName="CalSingleClusteringTool") ;

  const unsigned int nbDefaultCorrTools = 5 ;
  std::string defaultCorrTools[nbDefaultCorrTools] =
   { "LastLayerCorrTool", "ProfileTool", "CalValsCorrTool",
     "CalTkrLikelihoodTool", "CalTransvOffsetTool" } ;
  declareProperty("corrToolNames", m_corrToolNames
   = std::vector<std::string>(defaultCorrTools,defaultCorrTools+nbDefaultCorrTools)) ;

 }



StatusCode CalClustersAlg::initialize()

// This function does following initialization actions:
//    - extracts geometry constants from xml file using GlastDetSvc
//    - gets callNumber parameter from jobOptions file
//    - creates minimizer object (Midnight)
//    - sets pointer to the global function to be minimized
//      in Profile() function
//    - clears the global vector( g_elayer) with energies per layer in GeV

{
    if (CalReconActor::initialize(serviceLocator()).isFailure()) {
        return StatusCode::FAILURE ;
    }
    
    MsgStream log(msgSvc(), name()) ;
    StatusCode sc = StatusCode::SUCCESS ;
    
    sc = toolSvc()->retrieveTool(m_clusteringToolName,  m_clusteringTool);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to create " << m_clusteringToolName << endreq;
        return sc;
    }
    
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
    IEnergyCorr * tool ;
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

//StatusCode CalClustersAlg::retrieve()
//
//// Purpose:
////    - to get pointer to existing CalXtalRecCol
////    - to create new CalClusterCol and register it in TDS 
//
//{
//    StatusCode sc = StatusCode::SUCCESS;
//    MsgStream log(msgSvc(), name());  
//    
////    // get pointer to Calrecon directory in TDS
////    DataObject* pnode=0;
////    sc = eventSvc()->retrieveObject(EventModel::CalRecon::Event,pnode) ;
////
////    // if this directory doesn't yet exist - attempt to create it
////    if( sc.isFailure() ) {
////        sc = eventSvc()->registerObject(EventModel::CalRecon::Event,
////            new DataObject);
////                
////        // if can't create - put error message and return
////        if( sc.isFailure() ) {
////            log << MSG::ERROR << "Could not create CalRecon directory"
////                << endreq;
////            return sc;
////        }
////    }
//    
//    // attempt to get pointer to CalClusterCol, if it exists already
//    m_calClusterCol = SmartDataPtr<CalClusterCol> (eventSvc(),
//        EventModel::CalRecon::CalClusterCol);
//    
//    // if it doesn't exist - create it and register in TDS
//    if (!m_calClusterCol) {
//        m_calClusterCol = new CalClusterCol();
//        sc = eventSvc()->registerObject(EventModel::CalRecon::CalClusterCol,
//            m_calClusterCol);
//    }
//    
//    getKernel()->getClusters()->delClusters() ;
//    
//    return sc;
//}
//

//Purpose and method:
//
//   This function performs the calorimeter cluster reconstruction.
//   The main actions are:
//      - calculate energy sum
//                  energy per layer
//                  average position per layer
//                  quadratic spread per layer
//      - fit the particle direction using Fit_Direction() function
//      - calculate particle energy by profile fitting method
//          using Profile() function
//      - calculate particle energy by last layer correction method
//          using Leak() function
//      - store all calculated quantities in CalCluster object
// 
// TDS input: CalXtalRecCol
// TDS output: CalClustersCol

StatusCode CalClustersAlg::execute()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    log<<MSG::DEBUG<<"Begin"<<endreq ;
    getKernel()->reviewEvent() ;
    
    // non fatal errors
    // if there's no CalXtalRec then CalClustersAlg is not happening
  	if (!getKernel()->getXtalRecs())
      return StatusCode::SUCCESS;
    // other errors are fatal
    if (getKernel()->getStatus().isFailure())
      return StatusCode::FAILURE ;
      
    // call the Clustering tool
    log<<MSG::DEBUG<<"Clustering"<<endreq ;
    m_clusteringTool->findClusters() ;
    
    // loop over all found clusters
    log<<MSG::DEBUG<<"Corrections"<<endreq ;
    Event::CalClusterCol::const_iterator it ;
    for ( it = getKernel()->getClusters()->begin() ;
          it != getKernel()->getClusters()->end() ;
          ++it ) {       
        
        // if no tracker rec then fill slope from cluster
        if (getKernel()->getTkrNVertices()==0)
          getKernel()->setSlope((*it)->getDirection().z()) ;
 

        // apply corrections
        int itool = 0 ;
        std::vector<IEnergyCorr *>::const_iterator tool ;
        for ( tool = m_corrTools.begin() ;
              tool != m_corrTools.end() ;
              ++tool, ++itool ) {
            log<<MSG::DEBUG<<"Correction "<<itool<<endreq ;
            (*tool)->doEnergyCorr(*it);
        }
        
    }   // close loop on clusters
    
    // print the reconstruction results for debugging
    getKernel()->getClusters()->writeOut(log<<MSG::DEBUG);     
    
    log<<MSG::DEBUG<<"End"<<endreq ;
    return sc;
}

StatusCode CalClustersAlg::finalize()
 { return StatusCode::SUCCESS ; }




