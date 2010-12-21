
#include <CalRecon/ICalClassifyTool.h>
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

#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"
#include <algorithm>

/**   
* @class CalClassifyAlg
*
* @brief Algorithm for the classification of calorimeter clusters
*
* The idea is to use all the information available from the first pass
* of CalRecon, so containing CAL only information, in order to classify
* the clusters to tell what's their probability to be gamma/hadron/mip like,
* and potentially pass the most gamma-like cluster to the tracker reconstruction.
*
* Authors:  Johan Bregeon, Luca Baldini.
*/


bool SortByGamProb(Event::CalCluster *clusterA, Event::CalCluster *clusterB)
{
  return (clusterA->getGamProb() > clusterB->getGamProb());
}


class CalClassifyAlg : public Algorithm
{
public:
  
  //! Constructor
  CalClassifyAlg( const std::string & name, ISvcLocator * pSvcLocator ) ; 
  
  //! Destructor
  virtual ~CalClassifyAlg() {}; 
    
  virtual StatusCode initialize() ;
  
  StatusCode execute();
  
  StatusCode finalize() ;
        
private:
    
  /// Sort clusters according to the classification or not.
  bool              m_sortByClassification;
  /// name of Tool for finding clusters.
  StringProperty    m_classifierToolName;
  /// pointer to actual tool for finding clusters.
  ICalClassifyTool* m_classifierTool;
  //! package service.
  ICalReconSvc *    m_calReconSvc ;
} ;


//==============================================
// IMPLEMENTATION
//==============================================


DECLARE_ALGORITHM_FACTORY(CalClassifyAlg) ;

CalClassifyAlg::CalClassifyAlg(const std::string & name, ISvcLocator * pSvcLocator) :
  Algorithm(name,pSvcLocator),
  m_calReconSvc(0)
{   
  declareProperty("classifierToolName", m_classifierToolName = "CalClusterNBClassifyTool");
  declareProperty("sortByClassification", m_sortByClassification = true) ;
}

StatusCode CalClassifyAlg::initialize()
{
  MsgStream log(msgSvc(), name()) ;
  StatusCode sc = StatusCode::SUCCESS ;
  
  setProperties();
  
  if ( service("CalReconSvc",m_calReconSvc,true).isFailure() ) {
    log<<MSG::ERROR<<"Could not find CalReconSvc"<<endreq ;
    return StatusCode::FAILURE ;
  }

  sc = toolSvc()->retrieveTool(m_classifierToolName,  m_classifierTool);
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "  Unable to create " << m_classifierToolName << endreq;
    return sc;
  }
  
  return sc;
}

StatusCode CalClassifyAlg::execute()
{
  MsgStream log(msgSvc(), name());
  StatusCode sc = StatusCode::SUCCESS;
  
  log << MSG::DEBUG << "Begin CalClassifyAlg::execute()" << endreq ;

  // If there are no clusters, we should not classify them.
  // Check code below to get the cluster collectin through the SmartDataPtr -- TBD
  if ( !m_calReconSvc->getClusters() ) {
    return StatusCode::SUCCESS;
  }
  // Call the classifier tool!
  try {
    // Insure CalRecon/Event directory in TDS
    DataObject * pnode = 0 ;
    if ( (eventSvc()->retrieveObject(EventModel::CalRecon::Event,pnode)).isFailure()
	 && (eventSvc()->registerObject(EventModel::CalRecon::Event, new DataObject)).isFailure() ) { 
      throw CalException("Cannot register Event/CalRecon.") ;
    } 

    // Look up the CalClusterCol
    // TBD -- why do we call this SmartDataPtr stuff and not getClusters()
    Event::CalClusterCol* calClusterCol = 
      SmartDataPtr<Event::CalClusterCol>(eventSvc(),EventModel::CalRecon::CalClusterCol);

    // Call the tool to classify the clusters
    if ( m_classifierTool->classifyClusters(calClusterCol).isFailure() ) {
      sc = m_calReconSvc->handleError(name(),"classifier tool failure");
    }
       
    // Sort the cluster now that they are classified. Put the Uber cluster at the end though
    if( m_sortByClassification ) {
      log << MSG::DEBUG << "Sorting clusters using gam probability..." << endreq;
      sort (calClusterCol->begin(), calClusterCol->end()-1, SortByGamProb);
    }
    else {
      log << MSG::DEBUG << "Clusters are sorted by energy." << endreq;	  
    }
    
    // For debug
    log << MSG::DEBUG << "Last cluster is always the uber cluster:" << endreq;
    Event::CalClusterCol::const_iterator cluster;
    int clusterId = 1;
    for ( cluster = calClusterCol->begin(); cluster != calClusterCol->end(); cluster++)
      {
	log << MSG::DEBUG << "Full info for cluster (" << clusterId << ")\n" <<
	  **cluster << endreq;
	clusterId ++;
      }
    
  // Catch any exceptions here
  }
  catch( CalException & e ) {
    sc = m_calReconSvc->handleError(name()+" CalException",e.what());
  }
  catch( std::exception & e) {
    sc = m_calReconSvc->handleError(name()+" std::exception",e.what());
  }
  catch(...) {
    sc = m_calReconSvc->handleError(name(),"unknown exception");
  }
    
  log << MSG::DEBUG << "Done with CalClassifyAlg::execute()." << endreq;
  return sc;
}

StatusCode CalClassifyAlg::finalize()
{ 
  return StatusCode::SUCCESS; 
}
