
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/MsgStream.h"

#include "Event/TopLevel/EventModel.h"

#include "src/GCRRecon/IGcrReconTool.h"
#include "src/Utilities/GcrReconException.h"

//Needed by Engine4
#include "Event/TopLevel/Event.h"
//#include "Event/TopLevel/EventModel.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/GaudiException.h" 
#include "OnboardFilterTds/FilterStatus.h"
#include "OnboardFilterTds/ObfFilterStatus.h"
#include "OnboardFilterTds/FilterAlgTds.h"
#include "enums/TriggerBits.h"

#include "geometry/Point.h"
#include "geometry/Vector.h"

#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/GaudiException.h" 
#include "TkrUtil/ITkrGeometrySvc.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "Event/MonteCarlo/McParticle.h"


#include "Event/Recon/TkrRecon/TkrEventParams.h"



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

  StatusCode readGlastDet();  //read CAL geometry
  StatusCode getCalEntryExitPoints();  // if using TkrRecon, extrapolate Track1 initial dir to verify it enters the CAL

  ///gets the First McParticle launched for current event
  Event::McParticle* GcrReconAlg::findFirstMcParticle();


  //! Check if the event is a heavy ion candidate.
  bool isValidForGCR();

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
  
    //CAL dimensions parameters
  double             m_calZTop;
  double             m_calZBot;
  double             m_calXLo;
  double             m_calXHi;
  double             m_calYLo;
  double             m_calYHi;

  Point m_calEntryPoint;
  Point m_calExitPoint;

  //initial direction of the particle (initial Direction of TKR1 if using TKRTrack)
  Vector m_initDir;

    
} ;

namespace {
    int _noTkrColCount;
    bool _debug;
}

#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_ALGORITHM_FACTORY(GcrReconAlg) ;


GcrReconAlg::GcrReconAlg( const std::string & name, ISvcLocator * pSvcLocator ) : Algorithm(name,pSvcLocator), m_initDir(0.,0.,0.) 
{   
    // Declare the properties with these defaults
    declareProperty("GcrReconToolName", m_gcrReconToolName = "GcrReconTool");
    declareProperty("InitAxis", m_initAxis = "TKR");
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

    _noTkrColCount = 0;

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

    log<< MSG::DEBUG;
    _debug = (log.isActive());
    log << endreq;
    
    if (_debug) log<<MSG::DEBUG<<"GcrReconAlg::execute Begin"<<endreq ;
    
    // TEST:  Does TkrTrack enters CAL?
    //if (!(m_initAxis=="MC")){
      readGlastDet();   // Warning: code duplicated in GcrReconTool, should be passed as parameter of findGcrXtals
      
      if(getCalEntryExitPoints()==StatusCode::SUCCESS)   // Warning: code duplicated in GcrReconTool, should be passed as parameter of findGcrXtals
      {
	if((m_calEntryPoint.x()<m_calXLo) || (m_calEntryPoint.x()>m_calXHi) || (m_calEntryPoint.y()<m_calYLo) || (m_calEntryPoint.y()>m_calXHi)) 
	  { 
	    if(_debug) log<<MSG::DEBUG<<"track1 out of calorimeter"<<endreq ;

            return sc;
	  }
	else
	  if(_debug) log<<MSG::DEBUG<<"track1 in the calorimeter"<<endreq ; 

       }
       else{
         if(_debug) log<<MSG::DEBUG<<"no TKRtrack found"<<endreq ;
	 return sc;
       } 
      
      
    //}

    
    
    // Apply an event filter selection step, then find GCRs
    if(isValidForGCR())
      {
        Vector dir;
        Point pos;
        m_gcrReconTool->findGcrXtals(m_initAxis, m_calEntryPoint, m_calExitPoint, m_initDir);
      }

    if(_debug) log<<MSG::DEBUG<<"GcrReconAlg::execute End"<<endreq ;
    return sc;
}

StatusCode GcrReconAlg::finalize()
{ 
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG <<"finalize Begin" << endreq ;
    log << MSG::INFO  << _noTkrColCount << " event(s) without a TkrTrackCol " << endreq;
    log << MSG::DEBUG <<"finalize End"   << endreq ;
    return StatusCode::SUCCESS ; 
    
}

bool GcrReconAlg::isValidForGCR()
{
  MsgStream log(msgSvc(), name());
  if(m_HfcOrCnoTrig == "TriggerEng4"){
    log<<MSG::VERBOSE<<"@@@@@@@@ Using Trigger Engine 4"<<endreq ;
    SmartDataPtr<Event::EventHeader> pEvent(m_dataSvc, EventModel::EventHeader);
    unsigned int word2 =  ( pEvent==0? 0 : pEvent->triggerWordTwo());
    unsigned int Trig_gemengine = ((word2 >> enums::ENGINE_offset) & enums::ENGINE_mask);
    bool engine4ON2 = (Trig_gemengine==4);
    if(!engine4ON2){
      log<<MSG::INFO<<"@@@@@@@@ Trigger Engine 4 not set"<<endreq ;}
    bool passFilter=m_gcrReconTool->checkFilters(); // make a pass to compute the status word
    passFilter=false;//dummy line to avoid warning at compilation time
    return engine4ON2;
  }
  else{
    log<<MSG::VERBOSE<<"@@@@@@@@ Using Gamma,HFC,Mip,DGN OBF filters"<<endreq ;
    bool passFilter=m_gcrReconTool->checkFilters();
    
    return (passFilter);
 
  }
}

StatusCode GcrReconAlg::readGlastDet()
{
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());

  if(_debug) log << MSG::DEBUG << "GcrReconTool BEGIN readGlastDet()" << endreq ;  

  //TkrGeometrySvc used for access to tracker geometry info
  ITkrGeometrySvc*   m_geoSvc;
   // find TkrGeometrySvc service
  if (service("TkrGeometrySvc", m_geoSvc, true).isFailure()){
    log << MSG::ERROR << "Couldn't find the TkrGeometrySvc!" << endreq;
    return StatusCode::FAILURE;
  }
   
  m_calZTop = m_geoSvc->calZTop();
  m_calZBot = m_geoSvc->calZBot();
  
  double towerPitch = m_geoSvc->towerPitch();
  int xNum = m_geoSvc->numXTowers();
  int yNum = m_geoSvc->numYTowers();
  double calXWidth = m_geoSvc->calXWidth();
  double calYWidth = m_geoSvc->calYWidth();
  double deltaX = 0.5*(xNum*towerPitch - calXWidth);
  double deltaY = 0.5*(yNum*towerPitch - calYWidth);


  m_calXLo = m_geoSvc->getLATLimit(0,LOW)  + deltaX;
  m_calXHi = m_geoSvc->getLATLimit(0,HIGH) - deltaX;
  m_calYLo = m_geoSvc->getLATLimit(1,LOW)  + deltaY;
  m_calYHi = m_geoSvc->getLATLimit(1,HIGH) - deltaY;
  

  //  log << MSG::DEBUG << "Cal limits: YLo,YHi, XLo, XHi" << m_calYLo <<"," << m_calYHi << "," << m_calXLo << ","<< m_calXHi<< endreq;  
  if(_debug) log << MSG::DEBUG << "GcrReconTool END readGlastDet()" << endreq ;  
  return sc;
}



//---------------------------------------------------------------------------------

StatusCode GcrReconAlg::getCalEntryExitPoints(){
  //returns StatusCode::FAILURE if no TkrTrack found
  
  MsgStream log(msgSvc(), name());
  StatusCode sc = StatusCode::SUCCESS;
 
  Point initPos;
  Vector mcDir;


  if(m_initAxis == "MC"){
       if(_debug) log << MSG::DEBUG << " CAL Exit Entry points for initAxis==MC " << endreq;

       // ******the first McParticle initial position (when launched):
       Event::McParticle* firstMcParticle = findFirstMcParticle();

      const HepPoint3D& initPosition =firstMcParticle->initialPosition(); //initialPosition of firstMCParticle 
      // casting to Point of initPosition, necessary to define propagator steps

      if(firstMcParticle != 0)
	    mcDir = Vector(firstMcParticle->initialFourMomentum().vect().unit());
      else 
	    mcDir = Vector(-1e9,-1e9,-1e9);


      initPos.setX(initPosition.x());
      initPos.setY(initPosition.y());
      initPos.setZ(initPosition.z());


      // ******the first McParticle initial direction (when launched):
      m_initDir = Vector(mcDir.x(), mcDir.y(), mcDir.z());  
  }  
  
  else if(m_initAxis=="TKR"){
      if(_debug) log << MSG::DEBUG << " CAL Exit Entry points for initAxis==TKR " << endreq;

      //Locate and store a pointer to the data service
      DataSvc*           dataSvc;
      IService* iService = 0;
      if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
      throw GaudiException("Service [EventDataSvc] not found", name(), sc);
      dataSvc = dynamic_cast<DataSvc*>(iService);


      SmartDataPtr<Event::TkrTrackCol>   pTracks(dataSvc,EventModel::TkrRecon::TkrTrackCol);

	    if (pTracks){   
        	// Count number of tracks

        	int nTracks = pTracks->size();

        	if(_debug) log << MSG::DEBUG << "nTracks=" << nTracks << endreq;  
        	if(nTracks < 1) return StatusCode::FAILURE;
        	// Get the first Track - it should be the "Best Track"
        	Event::TkrTrackColConPtr pTrack = pTracks->begin();

		const Event::TkrTrack* track_1 = *pTrack;

		initPos = track_1->getInitialPosition();
		m_initDir = track_1->getInitialDirection().unit();

	    }
	    else
	{
        if(_noTkrColCount<5) {
            log << MSG::INFO << "no TkrTrackCol found : no subsequent processing" << endreq;
        }
        if(_noTkrColCount==4) {log << MSG::INFO << "further messages suppressed" << endreq;}
        _noTkrColCount++;
		return StatusCode::FAILURE;
	}
  }
  else
    {
      log<<MSG::ERROR<<"Invalid property "<<m_initAxis<<endreq;
      return StatusCode::FAILURE;
    }


  double x0 = initPos.x();
  double y0 = initPos.y();
  double z0 = initPos.z();

  double ux0 = m_initDir.x();
  double uy0 = m_initDir.y();
  double uz0 = m_initDir.z();


  
  if(_debug) log << MSG::DEBUG << "initPos=(" << x0 <<"," << y0 << "," << z0 <<")"<< endreq;  
  if(_debug) log << MSG::DEBUG << "m_initDir=(" << ux0 <<"," << uy0 << "," << uz0 <<")"<< endreq;  


  if (uz0!=0)
    {
      m_calEntryPoint = initPos+((m_calZTop-z0)/uz0)*m_initDir; 
      m_calExitPoint  = initPos+((m_calZBot-z0)/uz0)*m_initDir;        
      
    }
  else
    {
      if (ux0!=0)
        {
         
	  if (ux0>0){
	    m_calEntryPoint = initPos+((m_calXLo-x0)/ux0)*m_initDir; 
	    m_calExitPoint  = initPos+((m_calXHi-x0)/ux0)*m_initDir;
	  }
	  else{
	    m_calEntryPoint = initPos+((m_calXHi-x0)/ux0)*m_initDir; 
	    m_calExitPoint  = initPos+((m_calXLo-x0)/ux0)*m_initDir;
	    }
        }
      else if (uy0!=0)
        {
	  if (uy0>0){
	    m_calEntryPoint = initPos+((m_calYLo-y0)/uy0)*m_initDir; 
	    m_calExitPoint  = initPos+((m_calYHi-y0)/uy0)*m_initDir;
	    }
	  else{
	    m_calEntryPoint = initPos+((m_calYHi-y0)/uy0)*m_initDir; 
	    m_calExitPoint  = initPos+((m_calYLo-y0)/uy0)*m_initDir;
	    
	    }
        }
  
    }

  if(_debug) log << MSG::DEBUG << "CalEntryPoint=       " << '(' << m_calEntryPoint.x() << ',' << m_calEntryPoint.y() << ',' << m_calEntryPoint.z() << ')'  <<endreq;
  if(_debug) log << MSG::DEBUG << "CalExitPoint=       " << '(' << m_calExitPoint.x() << ',' << m_calExitPoint.y() << ',' << m_calExitPoint.z() << ')'  <<endreq;
    
  return sc;

}

//--------------------------------------------------------------------------

Event::McParticle* GcrReconAlg::findFirstMcParticle(){

  MsgStream log(msgSvc(), name());

  //m_log << MSG::INFO << "BEGIN findFirstMcParticle in GcrReconTool" << endreq;

  Event::McParticleCol* mcParticleCol = SmartDataPtr<Event::McParticleCol>(m_dataSvc, EventModel::MC::McParticleCol);
  // Event::McParticleCol::const_iterator mcParticleColIter=mcParticleCol->begin(); 

   
  Event::McParticle* firstMcParticle = 0; 
  if (mcParticleCol){  
    firstMcParticle = *mcParticleCol->begin();
    //m_log << MSG::INFO << "firstMcParticle particle ID= " << firstMcParticle->particleProperty()<< endreq;
  } 
  else 
    log << MSG::INFO << "no mcParticleCol" << endreq;
       
  //m_log << MSG::INFO << "firstMcParticle.x()= " << firstMcParticle->initialPosition().x() << endreq;

  //m_log << MSG::INFO << "END findFirstMcParticle in GcrReconTool" << endreq;
  
  return firstMcParticle;
   


}




