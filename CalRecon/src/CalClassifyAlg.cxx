
#include <CalRecon/ICalClassifyTool.h>
#include <CalRecon/ICalReconSvc.h>
#include <CalRecon/ICalClusteringTool.h>

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

#include "Clustering/StdClusterInfo.h"
#include "Clustering/MomentsClusterInfo.h"


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
  return (clusterA->getClassParams().getGamProb() > clusterB->getClassParams().getGamProb());
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
  bool              m_fullSortByClassification;
  float             m_maxEtoSortByClassification;
  /// name of Tool for finding clusters.
  StringProperty    m_classifierToolName;
  /// pointer to actual tool for finding clusters.
  ICalClassifyTool* m_classifierTool;
  // package service.
  ICalReconSvc *    m_calReconSvc;
  // Utility for filling clusters
  ICalClusterFiller* m_clusterInfo;
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
  declareProperty("maxEtoSortByClassification", m_maxEtoSortByClassification = 300.) ;
  declareProperty("fullSortByClassification", m_fullSortByClassification = false);
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

  // Cluster filling utility
  // -- To be replaced with a generic version soon
  m_clusterInfo = new MomentsClusterInfo(m_calReconSvc);

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

    // Get relation between clusters and xtals
    Event::CalClusterHitTabList* xTal2ClusTabList = SmartDataPtr<Event::CalClusterHitTabList>(eventSvc(),EventModel::CalRecon::CalClusterHitTab);
    Event::CalClusterHitTab* xTal2ClusTab = 0;
    if (xTal2ClusTabList) xTal2ClusTab = new Event::CalClusterHitTab(xTal2ClusTabList);
    Event::CalXtalRecCol* calXtalRecCol = SmartDataPtr<Event::CalXtalRecCol>(eventSvc(),EventModel::CalRecon::CalXtalRecCol);

    // Call the tool to classify the clusters
    if ( m_classifierTool->classifyClusters(calClusterCol).isFailure() ) {
      sc = m_calReconSvc->handleError(name(),"classifier tool failure");
    }

    // Cluster sorting (after classification). Put the Uber cluster at the end though
    // Only if energy of any clusters is smaller than some value.
    int   nCluEgtThr  = 0;
    float cluMaxGProb = 0.;
    int cluMaxGProbId = -1;
    for ( Event::CalClusterCol::const_iterator cluster = calClusterCol->begin();          
          cluster != calClusterCol->end()-1;          
          cluster++) {  
      if ((*cluster)->getMomParams().getEnergy()> m_maxEtoSortByClassification) {nCluEgtThr+=1;}
      if ((*cluster)->getClassParams().getGamProb()  > cluMaxGProb) {
        cluMaxGProb = (*cluster)->getClassParams().getGamProb();
        cluMaxGProbId = std::find(calClusterCol->begin(), calClusterCol->end(), *cluster) - calClusterCol->begin();   
      }
    }
    
    // iterator to be used in the rotate method below.
    Event::CalClusterCol::iterator cluMaxGProbIt ;
    cluMaxGProbIt = calClusterCol->begin();
    if (cluMaxGProbId > 0) std::advance(cluMaxGProbIt, cluMaxGProbId);

    if (nCluEgtThr==0){

      if (m_fullSortByClassification){
        log << MSG::DEBUG << "Clusters sorting: by gam probability." << endreq;
        sort (calClusterCol->begin(), calClusterCol->end()-1, SortByGamProb); }
      else { 
        log << MSG::DEBUG << "Clusters sorting: by energy and then putting the one with highest gam prob first." << endreq;
        rotate(calClusterCol->begin(), cluMaxGProbIt, cluMaxGProbIt+1 ); }
    }
    else {
      log << MSG::DEBUG << "Clusters sorting: by energy." << endreq;          
    }

    int numClusters = 0;
    if(calClusterCol) numClusters = calClusterCol->size();
    log << MSG::DEBUG << "There are " << numClusters << " clusters." << endreq;

    // new cluster uber2 = uber w/o 2nd cluster
    if(numClusters>1)
      {
      log << MSG::DEBUG << "More than one cluster : creating the Uber2 cluster." << endreq;
        // First get the cluster id for each xtal
        int myxtal2clus[16][8][12];
        int i,j,k;
        for(i=0;i<16;++i)
          for(j=0;j<8;++j)
            for(k=0;k<12;++k)
              myxtal2clus[i][j][k] = -1;
        
        Event::CalClusterCol::iterator clusIter = calClusterCol->begin();
        int iclu = 0;
        int nxtal2 = 0;
        while(clusIter != calClusterCol->end())
          {
            if(iclu>0 && iclu==numClusters-1) break;
            //
            Event::CalCluster* cl = *clusIter++;
            if(xTal2ClusTab)
              {
                std::vector<Event::CalClusterHitRel*> xTalRelVec = xTal2ClusTab->getRelBySecond(cl);
                //
                if(!xTalRelVec.empty())
                  {
                    std::vector<Event::CalClusterHitRel*>::const_iterator it = xTalRelVec.begin();
                    for (; it != xTalRelVec.end(); it++)
                      {
                        // get poiner to the reconstructed data for individual crystal
                        Event::CalXtalRecData* recData = (*it)->getFirst();
                        //
                        int itow=recData->getPackedId().getTower();
                        int ilay=recData->getPackedId().getLayer();
                        int icol=recData->getPackedId().getColumn();
                        myxtal2clus[itow][ilay][icol] = iclu;
                        if(iclu==1) ++nxtal2;
                      }
                  }
              }
            ++iclu;
          }
        log << MSG::DEBUG << "There are " << nxtal2 << " xtals in the second cluster." << endreq;
        //
        int nxtaluber2 = 0;
        // Get xtal list of xtals not in second cluster
        log << MSG::DEBUG << "Creating the list of xtals not in the second cluster." << endreq;
        Event::CalCluster* myuber2cluster = NULL;
        XtalDataList *xTalClus = new XtalDataList();
        if(calXtalRecCol)
          {
            for(Event::CalXtalRecCol::const_iterator xTalIter=calXtalRecCol->begin(); xTalIter != calXtalRecCol->end(); xTalIter++)
              {
                Event::CalXtalRecData* xTalData = *xTalIter;
                int itow=xTalData->getPackedId().getTower();
                int ilay=xTalData->getPackedId().getLayer();
                int icol=xTalData->getPackedId().getColumn();
                if(myxtal2clus[itow][ilay][icol]==1) continue;
                xTalClus->push_back(xTalData);
                ++nxtaluber2;
              }
          }
        log << MSG::DEBUG << "There are " << nxtaluber2 << " xtals not in the second cluster." << endreq;
        if(!xTalClus->empty())
          {
            // create and fill the cluster - from Tracy
            log << MSG::DEBUG << "Creating the Uber2 cluster with " << nxtaluber2 << " xtals." << endreq;
            myuber2cluster = m_clusterInfo->fillClusterInfo(xTalClus);
            
            std::string producerName("CalMSTClusteringTool/") ;
            producerName += myuber2cluster->getProducerName() ;
            myuber2cluster->setProducerName(producerName) ;
            myuber2cluster->clearStatusBit(Event::CalCluster::ALLXTALS);
            // Add cluster into the collection
            log << MSG::DEBUG << "Adding the Uber2 cluster into the collection." << endreq;
            calClusterCol->push_back(myuber2cluster);
            // reorder such that the last cluster is the uber cluster (and cluster uber2 is the next to last)
            log << MSG::DEBUG << "Reordering such that the last cluster is the uber cluster (and cluster uber2 is the next to last)." << endreq;
            rotate(calClusterCol->end()-2,calClusterCol->end()-1,calClusterCol->end());
          }
        if(xTalClus) delete xTalClus;
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

  if(xTal2ClusTab) delete xTal2ClusTab;

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
