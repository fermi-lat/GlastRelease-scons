
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/MsgStream.h"

#include "Event/TopLevel/EventModel.h"

#include "src/GCRSelect/IGcrSelectTool.h"
#include "src/Utilities/GcrSelectException.h"

/**   
* @class GcrSelectAlg
*
* @brief 
* 
*  
*/


class GcrSelectAlg : public Algorithm
{
public:
    //! constructor
    GcrSelectAlg( const std::string & name, ISvcLocator * pSvcLocator ); 
    
    //! destructor
    virtual ~GcrSelectAlg() {};
    
    virtual StatusCode initialize();

    StatusCode execute();

    StatusCode finalize() ;

private:

    
    //PRIVATE MEMBERS DATA
    //! correction tool names
    std::string      m_gcrSelectToolName ;
    
    //! correction tools
    IGcrSelectTool* m_gcrSelectTool ;
    
    //! package service
    //ICalSelectSvc *      m_calSelectSvc ;
} ;

#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_ALGORITHM_FACTORY(GcrSelectAlg) ;


GcrSelectAlg::GcrSelectAlg( const std::string & name, ISvcLocator * pSvcLocator ) : Algorithm(name,pSvcLocator)
{   
    // Declare the properties with these defaults
    declareProperty("GcrSelectToolName", m_gcrSelectToolName = "GcrSelectTool");
}



StatusCode GcrSelectAlg::initialize()
{
    // Purpose and Method: Initialize the algorithm:
    //                     - Initializes the vector of pointers to the various
    //                       energy correction tools to all for this particular
    //                       iteration of reconstruction

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Basic initialization first
    log << MSG::INFO << "GcrSelectAlg Initialization" << endreq;
    if( (sc = setProperties()).isFailure()) 
    {
        log << " GcrSelectAlg Initialization: didn't work!" << endreq;
        return sc;
    }
    log << endreq;
        
    

    if ((sc = toolSvc()->retrieveTool(m_gcrSelectToolName, m_gcrSelectTool)).isFailure())
    {
        log << MSG::ERROR << " Unable to create " << m_gcrSelectToolName << endreq ;
        return sc ;
    }

    return sc;
}


StatusCode GcrSelectAlg::execute()
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
    log<<MSG::INFO<<"GcrSelectAlg::execute Begin"<<endreq ;

   m_gcrSelectTool->selectGcrXtals();
  

    log<<MSG::INFO<<"GcrSelectAlg::execute End"<<endreq ;
    return sc;
}

StatusCode GcrSelectAlg::finalize()
{ 

    MsgStream log(msgSvc(), name());
    log<<MSG::INFO<<"GcrSelectAlg::finalize Begin"<<endreq ;
    log<<MSG::INFO<<"GcrSelectAlg::finalize End"<<endreq ;
   
    return StatusCode::SUCCESS ; 
    
}




