
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/MsgStream.h"

#include "Event/TopLevel/EventModel.h"

#include "src/GCRRecon/IGcrReconTool.h"
#include "src/Utilities/GcrReconException.h"

//#include <CalRecon/ICalReconSvc.h>

/**   
* @class GcrReconAlg
*
* @brief 
* 
*  
*/


class GcrReconAlg : public Algorithm
{
public:
    //! constructor
    GcrReconAlg( const std::string & name, ISvcLocator * pSvcLocator ); 
    
    //! destructor
    virtual ~GcrReconAlg() {};
    
    virtual StatusCode initialize();

    StatusCode execute();

    StatusCode finalize() ;

private:

    //! correction tool names
    std::string      m_gcrReconToolName ;
    
    //! correction tools
    IGcrReconTool* m_gcrReconTool ;
    
    //variable that indicates if we want to keep mcTrack direction or TrackReconTrack direction
    bool m_useMcDir;
} ;

#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_ALGORITHM_FACTORY(GcrReconAlg) ;


GcrReconAlg::GcrReconAlg( const std::string & name, ISvcLocator * pSvcLocator ) : Algorithm(name,pSvcLocator)
{   
    // Declare the properties with these defaults
    declareProperty("GcrReconToolName", m_gcrReconToolName = "GcrReconTool");
    declareProperty("UseMcDir", m_useMcDir = "true");
}



StatusCode GcrReconAlg::initialize()
{
    // Purpose and Method: Initialize the algorithm:
    //                     - Initializes the vector of pointers to the various
    //                       energy correction tools to all for this particular
    //                       iteration of reconstruction

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Basic initialization first
    log << MSG::INFO << "GcrReconAlg Initialization" << endreq;
    if( (sc = setProperties()).isFailure()) 
    {
        log << " GcrReconAlg Initialization: didn't work!" << endreq;
        return sc;
    }
    log << endreq;
        
    

    if ((sc = toolSvc()->retrieveTool(m_gcrReconToolName, m_gcrReconTool)).isFailure())
    {
        log << MSG::ERROR << " Unable to create " << m_gcrReconToolName << endreq ;
        return sc ;
    }

    return sc;
}


StatusCode GcrReconAlg::execute()
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
    
    log << MSG::INFO << "---------------@@@@@@@@@@@@@@ ------------" << endreq;
    log<<MSG::INFO<<"GcrReconAlg::execute Begin"<<endreq ;

    // Retrieve our TDS objects, we use Clusters to output corrected energy in CalEventEnergy
    
    // find GCRs
 
   m_gcrReconTool->findGcrXtals(m_useMcDir);
    /**try {
        if ((m_gcrReconTool->findGcrXtals()).isFailure()) {
            sc = m_calReconSvc->handleError(name(),"mip finding tool failure") ;
        }        
    } catch( CalException & e ) {
        sc = m_calReconSvc->handleError(name()+" CalException",e.what()) ;
    } catch( std::exception & e) {
        sc = m_calReconSvc->handleError(name()+" std::exception",e.what()) ;
    } catch(...) {
        sc = m_calReconSvc->handleError(name(),"unknown exception") ;
    }*/

    log<<MSG::INFO<<"GcrReconAlg::execute End"<<endreq ;
    return sc;
}

StatusCode GcrReconAlg::finalize()
{ 
    MsgStream log(msgSvc(), name());
    log<<MSG::INFO<<"GcrReconAlg::finalize Begin"<<endreq ;
    log<<MSG::INFO<<"GcrReconAlg::finalize End"<<endreq ;
    return StatusCode::SUCCESS ; 
    
}




