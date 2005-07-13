
#include <CalRecon/ICalClusteringTool.h>
#include <CalRecon/ICalReconSvc.h>

// for implementation
#include "src/Utilities/CalException.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

/**   
* @class CalClustersAlg
*
* @brief Algorithm for reconstruction of calo clusters
*
* The reconstruction here uses CalXtalRecCol to produce a CalClusterCol.
*  It evaluates the barycenter for each layer using the coordinates stored in
*  the CalXtalRecCol.
*
*/


class CalClustersAlg : public Algorithm
{
public:
    
    //! constructor
    CalClustersAlg( const std::string & name, ISvcLocator * pSvcLocator ) ; 
    //! destructor
    virtual ~CalClustersAlg() {}; 
    
    virtual StatusCode initialize() ;

        
    /*!Performs the reconstruction, creates one CalCluster object and stores
    * there the following results: 
    * - Energy per layer is computed and stored in CalCluster in MeV
    * - Barycenter per layer is also computed and stored in CalCluster
    * - The total energy, position and direction for whole cluster
    */        
    StatusCode execute();

    StatusCode finalize() ;
        
private:
    
    /// name of Tool for finding clusters
    StringProperty      m_clusteringToolName;

    /// pointer to actual tool for finding clusters
    ICalClusteringTool* m_clusteringTool;

    //! package service
    ICalReconSvc *      m_calReconSvc ;
} ;


//==============================================
// IMPLEMENTATION
//==============================================


DECLARE_ALGORITHM_FACTORY(CalClustersAlg) ;

CalClustersAlg::CalClustersAlg(const std::string & name, ISvcLocator * pSvcLocator )
                             : Algorithm(name,pSvcLocator), m_calReconSvc(0)
{   
    declareProperty ("clusteringToolName", m_clusteringToolName="CalSingleClusteringTool") ;
}

StatusCode CalClustersAlg::initialize()
{
    MsgStream log(msgSvc(), name()) ;
    StatusCode sc = StatusCode::SUCCESS ;
    
    if (service("CalReconSvc",m_calReconSvc,true).isFailure())
    {
        log<<MSG::ERROR<<"Could not find CalReconSvc"<<endreq ;
        return StatusCode::FAILURE ;
    }

    sc = toolSvc()->retrieveTool(m_clusteringToolName,  m_clusteringTool);
    if (sc.isFailure() ) 
    {
        log << MSG::ERROR << "  Unable to create " << m_clusteringToolName << endreq;
        return sc;
    }
        
    return sc;
}

StatusCode CalClustersAlg::execute()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    log<<MSG::DEBUG<<"Begin execute()"<<endreq ;
    
    // non fatal errors
    // if there's no CalXtalRec then CalClustersAlg is not happening
  	if (!m_calReconSvc->getXtalRecs())  return StatusCode::SUCCESS ;
      
    // call the clustering tool
    try {
        if (m_clusteringTool->findClusters(m_calReconSvc->getClusters()).isFailure()) {
            sc = m_calReconSvc->handleError(name(),"clustering tool failure") ;
        }        
    } catch( CalException & e ) {
        sc = m_calReconSvc->handleError(name()+" CalException",e.what()) ;
    } catch( std::exception & e) {
        sc = m_calReconSvc->handleError(name()+" std::exception",e.what()) ;
    } catch(...) {
        sc = m_calReconSvc->handleError(name(),"unknown exception") ;
    }
    
    log<<MSG::DEBUG<<"End execute()"<<endreq ;
    return sc;
}

StatusCode CalClustersAlg::finalize()
{ 
    return StatusCode::SUCCESS ; 
}




