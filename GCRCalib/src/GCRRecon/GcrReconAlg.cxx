
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
    
    //variable that indicates if we want to keep mcTrack direction or TrackReconTrack direction, default value: true
    bool m_useMcDir;
    //variable that indicates if we want to set HCF or CNOTrigger(CNO&LoCal&Tkr&RIO), default value: CnoTrig
    std::string m_HfcOrCnoTrig;
    
} ;

#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_ALGORITHM_FACTORY(GcrReconAlg) ;


GcrReconAlg::GcrReconAlg( const std::string & name, ISvcLocator * pSvcLocator ) : Algorithm(name,pSvcLocator)
{   
    // Declare the properties with these defaults
    declareProperty("GcrReconToolName", m_gcrReconToolName = "GcrReconTool");
    declareProperty("UseMcDir", m_useMcDir = "true");
    declareProperty("HFC_Or_TriggerEng4", m_HfcOrCnoTrig = "TriggerEng4");

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
    log << MSG::DEBUG << "GcrReconAlg Initialization" << endreq;
    if( (sc = setProperties()).isFailure()) 
    {
      log <<MSG::ERROR<< " GcrReconAlg Initialization: didn't work!" << endreq;
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
    
    log << MSG::DEBUG << "---------------@@@@@@@@@@@@@@ ------------" << endreq;
    log<<MSG::DEBUG<<"GcrReconAlg::execute Begin"<<endreq ;

    
    // Apply filter(TRIGGER Engine 4 "CNO Trigger" or HFC OnboardFilter), then find GCRs
       
   if(m_HfcOrCnoTrig == "TriggerEng4"){
      log<<MSG::DEBUG<<"@@@@@@@@ Using Trigger Engine 4"<<endreq ;
      if(m_gcrReconTool->TriggerEngine4ON()){
	 m_gcrReconTool->findGcrXtals(m_useMcDir);
      }
      else{
	 log<<MSG::DEBUG<<"@@@@@@@@ Trigger Engine 4 not set"<<endreq ;
      }    
   }
   else{
       log<<MSG::DEBUG<<"@@@@@@@@ Using HFC OBF"<<endreq ;
       if(!m_gcrReconTool->OBF_HFCVetoExist()){
	 m_gcrReconTool->findGcrXtals(m_useMcDir);
       }
       else{
	log<<MSG::DEBUG<<"@@@@@@@@ Vetoed Event"<<endreq ;
       }
   }

    log<<MSG::DEBUG<<"GcrReconAlg::execute End"<<endreq ;
    return sc;
}

StatusCode GcrReconAlg::finalize()
{ 
    MsgStream log(msgSvc(), name());
    log<<MSG::DEBUG<<"GcrReconAlg::finalize Begin"<<endreq ;
    log<<MSG::DEBUG<<"GcrReconAlg::finalize End"<<endreq ;
    return StatusCode::SUCCESS ; 
    
}




