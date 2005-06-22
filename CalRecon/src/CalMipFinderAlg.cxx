
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

#include "src/MipFinding/IMipFindingTool.h"

/**   
* @class CalMipFinderAlg
*
* @brief An algorithm for controlling and applying the various energy correction tools
*        used to determine the final event energy for GLAST
* 
* $Header$
*/


class CalMipFinderAlg : public Algorithm
{
public:
    //! constructor
    CalMipFinderAlg( const std::string & name, ISvcLocator * pSvcLocator ); 
    
    //! destructor
    virtual ~CalMipFinderAlg() {};
    
    virtual StatusCode initialize();

    StatusCode execute();

    StatusCode finalize() ;

private:

    //! correction tool names
    std::string      m_mipFinderName ;
    
    //! correction tools
    IMipFindingTool* m_mipFinder ;
} ;

#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_ALGORITHM_FACTORY(CalMipFinderAlg) ;


CalMipFinderAlg::CalMipFinderAlg( const std::string & name, ISvcLocator * pSvcLocator )
 : Algorithm(name,pSvcLocator)
{   
    // Declare the properties with these defaults
    declareProperty("MipFinderName", m_mipFinderName = "StdMipFindingTool");
}



StatusCode CalMipFinderAlg::initialize()
{
    // Purpose and Method: Initialize the algorithm:
    //                     - Initializes the vector of pointers to the various
    //                       energy correction tools to all for this particular
    //                       iteration of reconstruction

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Basic initialization first
    log << MSG::INFO << "CalMipFinderAlg Initialization";
    if( (sc = setProperties()).isFailure()) 
    {
        log << " didn't work!" << endreq;
        return sc;
    }
    log << endreq;
        
    if ((sc = toolSvc()->retrieveTool(m_mipFinderName, m_mipFinder)).isFailure())
    {
        log << MSG::ERROR << " Unable to create " << m_mipFinderName << endreq ;
        return sc ;
    }

    return sc;
}


StatusCode CalMipFinderAlg::execute()
{
    //Purpose and method: Primary driver for running the various energy correction tools
    //                    Also creates and registers the output TDS object for the Event Energy
    //                    and retrieves the "best" energy as the one to use this event.
    // 
    // TDS input:  CalClusterCol
    // TDS output: CalEventEnergy (with a vector of CalCorToolResult objects for each tool run)
    //
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    log<<MSG::DEBUG<<"Begin"<<endreq ;

    // Retrieve our TDS objects, we use Clusters to output corrected energy in CalEventEnergy
    //Event::CalClusterCol*  calClusters = SmartDataPtr<Event::CalClusterCol>(eventSvc(),EventModel::CalRecon::CalClusterCol);
    //Event::CalEventEnergy* calEnergy   = SmartDataPtr<Event::CalEventEnergy>(eventSvc(),EventModel::CalRecon::CalEventEnergy);
    //Event::TkrVertexCol*   tkrVertices = SmartDataPtr<Event::TkrVertexCol>(eventSvc(),EventModel::TkrRecon::TkrVertexCol);

    // find mips
    sc = m_mipFinder->findMIPCandidates();

    log<<MSG::DEBUG<<"End"<<endreq ;
    return sc;
}

StatusCode CalMipFinderAlg::finalize()
{ 
    return StatusCode::SUCCESS ; 
}




