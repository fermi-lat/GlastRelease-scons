/** @file CalClustersAlg.cxx
    @brief Implementation of CalClustersAlg

*/
#include "CalClustersAlg.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

#include "EnergyCorr.h"   // for special downcast

static const AlgFactory<CalClustersAlg>  Factory;
const IAlgFactory& CalClustersAlgFactory = Factory;

using namespace Event;

CalClustersAlg::CalClustersAlg(const std::string& name,
                               ISvcLocator* pSvcLocator):
Algorithm(name, pSvcLocator)
{
    
    // declaration of parameter needed to distinguish 2 calls
    // of CalClustersAlg:
    // 1st - before TkrRecon and 2nd - after TkrRecon
    
    declareProperty("callNumber",m_callNumber=0);
    declareProperty ("clusterToolName", m_clusterToolName="SingleClusterTool");
    declareProperty ("lastLayerToolName", m_lastLayerToolName="LastLayerCorrTool");
    declareProperty ("profileToolName", m_profileToolName="ProfileTool");
    declareProperty ("calValsCorrToolName", m_calValsCorrToolName="CalValsCorrTool");
    
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
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    
    
    // get pointer to GlastDetSvc
    sc = service("GlastDetSvc", detSvc);
    
    // if GlastDetSvc isn't available - put error message and return
    if(sc.isFailure())
    {
        log << MSG::ERROR << "GlastDetSvc could not be found" <<endreq;
        return sc;
    }
    
    
    // extracting detector geometry constants from xml file
    
    double value;
    if(!detSvc->getNumericConstByName(std::string("CALnLayer"), &value)) 
    {
        log << MSG::ERROR << " constant " << " CALnLayer "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } else m_CalnLayers = int(value);
    
    if(!detSvc->getNumericConstByName(std::string("CsIWidth"),&m_CsIWidth))
    {
        log << MSG::ERROR << " constant " << " CsIWidth "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } 
    
    if(!detSvc->getNumericConstByName(std::string("CsIHeight"),&m_CsIHeight))
    {
        log << MSG::ERROR << " constant " << " CsIHeight "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } 
    
    
    
    
    // get callNumber parameter from jobOptions file    
    setProperties();
    log << MSG::INFO << "CalClustersAlg: callNumber = " 
        << m_callNumber << endreq;
    
    // set global constants with geometry parameters (in cm)
    // used by Profile() function
    
    sc = toolSvc()->retrieveTool(m_clusterToolName,  m_clusterTool);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to create " << m_clusterToolName << endreq;
        return sc;
    }
    
    sc = toolSvc()->retrieveTool(m_lastLayerToolName, m_lastLayerTool);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to create " << m_lastLayerToolName << endreq;
        return sc;
    }

    sc = toolSvc()->retrieveTool(m_profileToolName, m_profileTool);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to create " << m_profileToolName << endreq;
        return sc;
    }

    sc = toolSvc()->retrieveTool(m_calValsCorrToolName,m_calValsCorrTool);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to create " << m_calValsCorrToolName << endreq;
        return sc;
    }
    return sc;
}

StatusCode CalClustersAlg::retrieve()

// Purpose:
//    - to get pointer to existing CalXtalRecCol
//    - to create new CalClusterCol and register it in TDS 
//
//  TDS input:   CalXtalRecCol
//  TDS ourput  CalClusterCol

{
    
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    
    DataObject* pnode=0;
    
    
    // get pointer to Calrecon directory in TDS
    sc = eventSvc()->retrieveObject(EventModel::CalRecon::Event, pnode );
    
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
    if (!m_calClusterCol )
    {
        m_calClusterCol = 0;
        m_calClusterCol = new CalClusterCol();
        sc = eventSvc()->registerObject(EventModel::CalRecon::CalClusterCol,
            m_calClusterCol);
    }
    
    // if it exists - delete all clusters
    else
    {
        m_calClusterCol->delClusters();
    }
    
    m_clusterTool->setClusterCol(m_calClusterCol);
    
    // get pointer to CalXtalRecCol
    m_calXtalRecCol = SmartDataPtr<CalXtalRecCol>(eventSvc(),
        EventModel::CalRecon::CalXtalRecCol); 
    
    
    return sc;
}

StatusCode CalClustersAlg::execute()

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


{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    //get pointers to the TDS data structures
    sc = retrieve();
    
//    const Point p0(0.,0.,0.);
    
    // variable indicating ( if >0) the presence of tracker
    // reconstruction output
    int rectkr=0;  
    
    int ntracks=0;
    Vector trackDirection;
    Point trackVertex;
    double slope;
    
    
    // get pointer to the tracker vertex collection
    SmartDataPtr<TkrVertexCol> tkrRecData(eventSvc(),
        EventModel::TkrRecon::TkrVertexCol);
    
    // if reconstructed tracker data doesn't exist - put the debugging message
    if (tkrRecData == 0) {
        log << MSG::DEBUG << "No TKR Reconstruction available " << endreq;
        // return sc;
    }
    
    
    // if data exist and number of tracks not zero 
    // - get information of track position and direction 
    else
    {
        // First get reconstructed direction from tracker
        ntracks = tkrRecData->size();
        log << MSG::DEBUG << "number of tracks = " << ntracks << endreq;
        
        
        if (ntracks > 0) {
            rectkr++;
            trackDirection = tkrRecData->front()->getDirection();
            trackVertex = tkrRecData->front()->getPosition();
            slope = fabs(trackDirection.z());
            log << MSG::DEBUG << "track direction = " << slope << endreq;
            
        } else {
            log << MSG::DEBUG << "No reconstructed tracks " << endreq;
        }	
    }
    
    
    // call the Cluster tool to find clusters
    m_clusterTool->findClusters(m_calXtalRecCol);
    
    
    // loop over all found clusters
    for (Event::CalClusterCol::const_iterator it = m_calClusterCol->begin();
    it != m_calClusterCol->end(); it++){       
        
        // if no tracker rec then fill slope from cluster
        if(!rectkr) slope = (*it)->getDirection().z();
 
        // do last layer correlation method

        m_lastLayerTool->setTrackSlope(slope);
        m_lastLayerTool->doEnergyCorr((*it)->getEnergySum(),(*it));

        // eleak is observed + estimated leakage energy
        double eleak = m_lastLayerTool->getEnergyCorr() + (*it)->getEnergySum();

        // iteration
        m_lastLayerTool->doEnergyCorr(eleak,(*it));       
        eleak = m_lastLayerTool->getEnergyCorr() + (*it)->getEnergySum();

        
        // Do profile fitting - use StaticSlope because of static functions
        // passed to minuit fitter

        dynamic_cast<EnergyCorr*>(m_profileTool)->setStaticSlope(slope);
        m_profileTool->doEnergyCorr((*it)->getEnergySum(),(*it));
        
         // get corrections from CalValsTool... self contained

        m_calValsCorrTool->doEnergyCorr((*it)->getEnergySum(),(*it));

        // calculating the transverse offset of average position in the calorimeter
        // with respect to the position predicted from tracker information
        double calTransvOffset = 0.;
        if(ntracks>0){
            Vector calOffset = ((*it)->getPosition()) - trackVertex;
            double calLongOffset = trackDirection*calOffset;
            calTransvOffset =sqrt(calOffset.mag2() - calLongOffset*calLongOffset);
            
        }
        
        // store the calculated quantities back in this CalCluster object. Note
        // temporary ugly kluge of overwriting all but two existing elements!
        (*it)->initialize(eleak,
            (*it)->getEneLayer(),
            (*it)->getPosLayer(),
            (*it)->getRmsLayer(),
            (*it)->getRmsLong(),
            (*it)->getRmsTrans(),
            (*it)->getDirection(),
            calTransvOffset);
    }   // close loop on clusters
    
    // print the reconstruction results for debugging
    m_calClusterCol->writeOut(log << MSG::DEBUG);
    
    return sc;
}

StatusCode CalClustersAlg::finalize()
{
    StatusCode sc = StatusCode::SUCCESS;
        
    return sc;
}




