
// Tool and Gaudi related stuff
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/DeclareFactoryEntries.h"

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
    ICalReconSvc*     m_calReconSvc;

    //! Utility for filling clusters
    ICalClusterFiller*   m_clusterInfo;
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
    XtalDataVec* xTalClus = new XtalDataVec();

    xTalClus->clear();
    Event::CalXtalRecCol::const_iterator it ;
    for (Event::CalXtalRecCol::const_iterator it = m_calReconSvc->getXtalRecs()->begin() ; 
                it != m_calReconSvc->getXtalRecs()->end(); ++it )
    {
        // get pointer to the reconstructed data for given crystal
	    Event::CalXtalRecData * recData = *it ;
        xTalClus->push_back(recData) ;
    }

    Event::CalCluster* cluster = m_clusterInfo->fillClusterInfo(xTalClus);
	cluster->setStatusBit(Event::CalCluster::ALLXTALS); 

    calClusterCol->push_back(cluster);

    delete xTalClus;
  
    return StatusCode::SUCCESS ;
}

