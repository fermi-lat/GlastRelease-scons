
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
  declareProperty("callNumber",m_callNumber=0) ;
    
  declareProperty ("clusteringToolName", m_clusteringToolName="CalSingleClusteringTool") ;
  declareProperty ("lastLayerToolName", m_lastLayerToolName="LastLayerCorrTool") ;
  declareProperty ("profileToolName", m_profileToolName="ProfileTool") ;
  declareProperty ("calValsCorrToolName", m_calValsCorrToolName="CalValsCorrTool") ;
  // disabled for now
  //declareProperty ("calTkrLikelihoodToolName", m_tkrLikelihoodToolName="CalTkrLikelihoodTool") ;
  declareProperty ("calTkrLikelihoodToolName", m_tkrLikelihoodToolName="") ;
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
    MsgStream log(msgSvc(), name()) ;
    StatusCode sc = StatusCode::SUCCESS ;
    
    m_data = new CalClusteringData(serviceLocator()) ;
    if (m_data->getStatus()==StatusCode::FAILURE)
     { return StatusCode::FAILURE ; }
        
    // get callNumber parameter from jobOptions file    
    setProperties();
    
    log << MSG::INFO << "CalClustersAlg: callNumber = " 
        << m_callNumber << endreq;
    
    // set global constants with geometry parameters (in cm)
    // used by Profile() function
    
    sc = toolSvc()->retrieveTool(m_clusteringToolName,  m_clusteringTool);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to create " << m_clusteringToolName << endreq;
        return sc;
    }
        
    if ( m_lastLayerToolName.value() != "" ) {
        sc = toolSvc()->retrieveTool(m_lastLayerToolName, m_lastLayerTool);
        if (sc.isFailure() ) {
            log << MSG::ERROR << "  Unable to create " << m_lastLayerToolName << endreq;
            return sc;
        }
    }
    else {
        m_lastLayerTool = 0 ;
    }

    if ( m_profileToolName.value() != "" ) {
        sc = toolSvc()->retrieveTool(m_profileToolName, m_profileTool);
        if (sc.isFailure() ) {
            log << MSG::ERROR << "  Unable to create " << m_profileToolName << endreq;
            return sc;
        }
    }
    else {
        m_profileTool = 0 ;
    }

    if ( m_calValsCorrToolName.value() != "" ) {
        sc = toolSvc()->retrieveTool(m_calValsCorrToolName,m_calValsCorrTool);
        if (sc.isFailure() ) {
            log << MSG::ERROR << "  Unable to create " << m_calValsCorrToolName << endreq;
            return sc;
        }
    }
    else {
        m_calValsCorrTool = 0 ;
    }

    if ( m_tkrLikelihoodToolName.value() != "" ) {
      printf("\033[32;41;1mhello \033[0m\n");
        sc = toolSvc()->retrieveTool(m_tkrLikelihoodToolName,m_tkrLikelihoodTool);
        if (sc.isFailure() ) {
            log << MSG::ERROR << "  Unable to create " << m_tkrLikelihoodToolName << endreq;
            return sc;
        }
    }
    else {
        m_calValsCorrTool = 0 ;
    }
    return sc;
}

StatusCode CalClustersAlg::retrieve()

// Purpose:
//    - to get pointer to existing CalXtalRecCol
//    - to create new CalClusterCol and register it in TDS 

{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());  
    
    // get pointer to Calrecon directory in TDS
    DataObject* pnode=0;
    sc = eventSvc()->retrieveObject(EventModel::CalRecon::Event,pnode) ;

    // if this directory doesn't yet exist - attempt to create it
    if( sc.isFailure() ) {
        sc = eventSvc()->registerObject(EventModel::CalRecon::Event,
            new DataObject);
                
        // if can't create - put error message and return
        if( sc.isFailure() ) {
            log << MSG::ERROR << "Could not create CalRecon directory"
                << endreq;
            return sc;
        }
    }
    
    // attempt to get pointer to CalClusterCol, if it exists already
    m_calClusterCol = SmartDataPtr<CalClusterCol> (eventSvc(),
        EventModel::CalRecon::CalClusterCol);
    
    // if it doesn't exist - create it and register in TDS
    if (!m_calClusterCol) {
        m_calClusterCol = new CalClusterCol();
        sc = eventSvc()->registerObject(EventModel::CalRecon::CalClusterCol,
            m_calClusterCol);
    }
    
    // if it exists - delete all clusters
    else {
        m_calClusterCol->delClusters();
    }
    
    // get pointer to CalXtalRecCol
    m_calXtalRecCol = SmartDataPtr<CalXtalRecCol>(eventSvc(),
        EventModel::CalRecon::CalXtalRecCol); 
    if (!m_calXtalRecCol) {
        log<<MSG::VERBOSE<<"No CalXtalRecCol"<<endreq ;
    }
        
    return sc;
}


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
    
    // get pointers to the TDS data structures
    sc = retrieve() ;
    // non fatal errors
    // if there's no CalXtalRec then CalClustersAlg is not happening
  	if (!m_calXtalRecCol)
      return StatusCode::SUCCESS;
    // other retreive errors are fatal
    if (sc.isFailure())
      return StatusCode::FAILURE;
         
    // update otherinputs
    m_data->beginEvent() ;
            
    // call the Clustering tool to find clusters
    m_clusteringTool->findClusters(m_calXtalRecCol,m_calClusterCol);
    
    // loop over all found clusters
    Event::CalClusterCol::const_iterator it ;
    for ( it = m_calClusterCol->begin() ;
          it != m_calClusterCol->end() ;
          ++it ) {       
        
        // if no tracker rec then fill slope from cluster
        if (m_data->getTkrNVertices()==0)
          m_data->setSlope((*it)->getDirection().z()) ;
 
        // do last layer correlation method
        double elastlayer = 0 ;
        if (m_lastLayerTool) {
            
            m_lastLayerTool->doEnergyCorr(m_data,*it) ;
            
            // eleak is observed + estimated leakage energy
            // double elastlayer = m_lastLayerTool->getEnergyCorr() + (*it)->getEnergySum();
            elastlayer = m_lastLayerTool->getEnergyCorr();

            // iteration commented out, now done in LastLayerTool.cxx
            //m_lastLayerTool->doEnergyCorr(elastlayer,(*it));       
            //elastlayer = m_lastLayerTool->getEnergyCorr() + (*it)->getEnergySum() ;
        }
        
        
        // [Phillipe&?] Do profile fitting - use StaticSlope because of static functions
        // passed to minuit fitter
        if (m_profileTool) {
//            dynamic_cast<EnergyCorr*>(m_profileTool)->setStaticSlope(slope);
//            m_profileTool->setTrackSlope(slope);
            m_profileTool->doEnergyCorr(m_data,*it);
        }
        
        // [Bill Atwood] get corrections from CalValsTool... self contained
        if (m_tkrLikelihoodTool) {       
            m_tkrLikelihoodTool->doEnergyCorr(m_data,*it);
        }

        // [Pol] get corrections from CalValsTool... self contained
        if (m_calValsCorrTool) {       
            m_calValsCorrTool->doEnergyCorr(m_data,*it);
        }

        // calculating the transverse offset of average position in the calorimeter
        // with respect to the position predicted from tracker information
        double calTransvOffset = 0.;
        if (m_data->getTkrNVertices()>0) {
            Vector calOffset = ((*it)->getPosition()) - m_data->getTkrFrontVertexPosition() ;
            double calLongOffset = m_data->getTkrFrontVertexDirection()*calOffset;
            calTransvOffset = sqrt(calOffset.mag2() - calLongOffset*calLongOffset);
        }
        
        // store the calculated quantities back in this CalCluster object. Note
        // temporary ugly kluge of overwriting all but two existing elements!
        (*it)->initialize(elastlayer, // usefull ?
            (*it)->getEneLayer(),
            (*it)->getPosLayer(),
            (*it)->getRmsLayer(),
            (*it)->getRmsLong(),
            (*it)->getRmsTrans(),
            (*it)->getDirection(),
            calTransvOffset) ;
    }   // close loop on clusters
    
    // print the reconstruction results for debugging
    m_calClusterCol->writeOut(log<<MSG::DEBUG);     
    
    return sc;
}

StatusCode CalClustersAlg::finalize()
 { return StatusCode::SUCCESS ; }




