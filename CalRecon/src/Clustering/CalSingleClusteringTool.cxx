
// Tool and Gaudi related stuff
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/DeclareFactoryEntries.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"

#include <CalRecon/ICalReconSvc.h>
#include <CalRecon/ICalClusteringTool.h>
#include "StdClusterInfo.h"
#include "MomentsClusterInfo.h"

class CalSingleClusteringTool : public AlgTool, virtual public ICalClusteringTool
{
public :
  
    /// Standard Gaudi Tool interface constructor
    CalSingleClusteringTool(const std::string& type,
                            const std::string& name,
                            const IInterface* parent );

    virtual ~CalSingleClusteringTool() {};
    
    /// @brief Intialization of the tool
    virtual StatusCode initialize() ;

    /// @brief Default cluster finding framework
    virtual StatusCode findClusters(Event::CalClusterCol* calClusterCol) ;

private:
    //! Service for basic Cal info
    ICalReconSvc*      m_calReconSvc;

    //! Event Service member directly useable by concrete classes.
    IDataProviderSvc*  m_dataSvc;

    //! Utility for filling clusters
    ICalClusterFiller* m_clusterInfo;
} ;

DECLARE_TOOL_FACTORY(CalSingleClusteringTool) ;

CalSingleClusteringTool::CalSingleClusteringTool(const std::string & type, 
                                                 const std::string & name, 
                                                 const IInterface* parent)
                                               : AlgTool(type,name,parent)
{ 
    declareInterface<ICalClusteringTool>(this) ; 
}
    
StatusCode CalSingleClusteringTool::initialize()
{
    MsgStream log(msgSvc(),name()) ;
    StatusCode sc = StatusCode::SUCCESS;

    if ((sc = service("CalReconSvc",m_calReconSvc,true)).isFailure())
    {
        throw GaudiException("Service [CalReconSvc] not found", name(), sc);
    }

    if(service( "EventDataSvc", m_dataSvc, true ).isFailure()) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    // Cluster filling utility
    // -- To be replaced with a generic version soon
    m_clusterInfo = new MomentsClusterInfo(m_calReconSvc);

    return StatusCode::SUCCESS ;
}

StatusCode CalSingleClusteringTool::findClusters(Event::CalClusterCol* calClusterCol)
{
    //Purpose and method:
    //
    //   This function performs the calorimeter cluster reconstruction.
    //   The main actions are:
    //      - calculate energy sum
    //                  energy per layer
    //                  average position per layer
    //                  quadratic spread per layer
    //      - fit the particle direction using Fit_Direction() function
    //      - store all calculated quantities in CalCluster objects
    // 
    // TDS input: CalXtalRecCol
    // TDS output: CalClustersCol
    // prepare the initital set of xtals
    XtalDataList* xTalClus = new XtalDataList();

    // Create a Xtal to Cluster relations list
    Event::CalClusterHitTabList* xTal2ClusTabList = new Event::CalClusterHitTabList();
    xTal2ClusTabList->clear();

    // Register the list in the TDS (which will assume ownership of the list, but not the table)
    if (m_dataSvc->registerObject(EventModel::CalRecon::CalClusterHitTab, xTal2ClusTabList).isFailure())
    {
        throw GaudiException("Unable to register xTal to Cluster table in TDS", name(), StatusCode::FAILURE);
    }

    xTalClus->clear();
    for (Event::CalXtalRecCol::const_iterator it = m_calReconSvc->getXtalRecs()->begin() ; 
                it != m_calReconSvc->getXtalRecs()->end(); ++it )
    {
        // get pointer to the reconstructed data for given crystal
        Event::CalXtalRecData * recData = *it ;
        xTalClus->push_back(recData) ;
    }

    calClusterCol->clear() ;

    // Get the cluster instance
    Event::CalCluster* cluster = m_clusterInfo->fillClusterInfo(xTalClus);

    std::string producerName("CalSingleClusteringTool/") ;
    producerName += cluster->getProducerName() ;
    cluster->setProducerName(producerName) ;
    cluster->setStatusBit(Event::CalCluster::ALLXTALS); 
    calClusterCol->push_back(cluster);

    // Loop through again to make the relations
    for (Event::CalXtalRecCol::const_iterator it = m_calReconSvc->getXtalRecs()->begin() ; 
                it != m_calReconSvc->getXtalRecs()->end(); ++it )
    {
        // get pointer to the reconstructed data for given crystal
        Event::CalXtalRecData * recData = *it ;

        // Even though only one cluster, we still need to make the relations!
        Event::CalClusterHitRel* xTal2ClusRel = new Event::CalClusterHitRel(recData,cluster);

        xTal2ClusTabList->push_back(xTal2ClusRel);
    }

    delete xTalClus;
  
    return StatusCode::SUCCESS ;
}

