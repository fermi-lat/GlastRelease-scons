
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/MsgStream.h"

#include "Event/TopLevel/EventModel.h"

#include "src/GCRRecon/IGcrReconTool.h"
#include "src/Utilities/GcrReconException.h"

//Needed by Engine4
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/GaudiException.h" 
#include "OnboardFilterTds/FilterStatus.h"
#include "OnboardFilterTds/ObfFilterStatus.h"
#include "OnboardFilterTds/FilterAlgTds.h"
#include "enums/TriggerBits.h"

#include "geometry/Point.h"
#include "geometry/Vector.h"

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
  //! Check if the event is a heavy ion candidate.
  bool passTrigger();
  void selectEventAxis(Vector&,Point&);

    //! correction tool names
    std::string      m_gcrReconToolName ;
    
    //! correction tools
    IGcrReconTool* m_gcrReconTool ;
    
  ///variable that indicates if we want to keep mcTrack 
  ///direction or TrackReconTrack direction, default value: true
    std::string m_initAxis;
  ///variable that indicates if we want to set HCF or 
  ///CNOTrigger(CNO&LoCal&Tkr&RIO), default value: CnoTrig
    std::string m_HfcOrCnoTrig;

  /// Pointer to the Gaudi data provider service
  DataSvc*           m_dataSvc;
    
} ;

#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_ALGORITHM_FACTORY(GcrReconAlg) ;


GcrReconAlg::GcrReconAlg( const std::string & name, ISvcLocator * pSvcLocator ) : Algorithm(name,pSvcLocator)
{   
    // Declare the properties with these defaults
    declareProperty("GcrReconToolName", m_gcrReconToolName = "GcrReconTool");
    declareProperty("InitAxis", m_initAxis = "MC");
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
        
  //Locate and store a pointer to the data service
  IService* iService = 0;
  if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
    throw GaudiException("Service [EventDataSvc] not found", name(), sc);
  m_dataSvc = dynamic_cast<DataSvc*>(iService);
    

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
    if(passTrigger())
      {
        Vector dir;
        Point pos;
        selectEventAxis(dir, pos);
        //m_gcrReconTool->findGcrXtals(m_initAxis, dir, pos);
        m_gcrReconTool->findGcrXtals(m_initAxis);
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

bool GcrReconAlg::passTrigger()
{
  MsgStream log(msgSvc(), name());
  if(m_HfcOrCnoTrig == "TriggerEng4"){
    log<<MSG::DEBUG<<"@@@@@@@@ Using Trigger Engine 4"<<endreq ;
    SmartDataPtr<Event::EventHeader> pEvent(m_dataSvc, EventModel::EventHeader);
    unsigned int word2 =  ( pEvent==0? 0 : pEvent->triggerWordTwo());
    unsigned int Trig_gemengine = ((word2 >> enums::ENGINE_offset) & enums::ENGINE_mask);
    bool engine4ON2 = (Trig_gemengine==4);
    if(!engine4ON2){
      log<<MSG::DEBUG<<"@@@@@@@@ Trigger Engine 4 not set"<<endreq ;}
    return engine4ON2;
    //    return m_gcrReconTool->TriggerEngine4ON();
  }
  else{
    log<<MSG::DEBUG<<"@@@@@@@@ Using HFC OBF"<<endreq ;
    bool vetoExists=true;
    
    SmartDataPtr<OnboardFilterTds::ObfFilterStatus> 
      obfStatus(m_dataSvc, "/Event/Filter/ObfFilterStatus");
    if (obfStatus)
      {
        // Pointer to our retrieved objects
        const OnboardFilterTds::IObfStatus* obfResult = 0;
        obfResult = obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::HFCFilter);
        if(obfResult){
          unsigned int statusHFC32 = obfResult->getStatus32();
          unsigned int vetoHFC= obfResult->getVetoBit();
          vetoExists = (statusHFC32 & vetoHFC)>0;   
          log << MSG::INFO << "(statusHFC32 & vetoHFC)>0= " << vetoExists << endreq;
        }     else{
          log << MSG::INFO <<  "no obfResult" << endreq;}
      }
    else{
      log << MSG::INFO << "no obfStatus"<< endreq;}
    
    return (!vetoExists);
    
    //	log<<MSG::DEBUG<<"@@@@@@@@ Vetoed Event"<<endreq ;
  }
}

void GcrReconAlg::selectEventAxis(Vector &dir, Point &pos)
{
  return;
}


