// $Header$
// 
//  Original author: Toby Burnett tburnett@u.washington.edu
//

#include "GuiSvc/GuiSvc.h"

#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"

#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/Property.h"
#include "GaudiKernel/IAppMgrUI.h"
#include "GaudiKernel/SmartIF.h"
#include "GaudiKernel/IObjManager.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/IToolFactory.h"

#include "GuiSvc/IGuiTool.h"

// includes for implementing GuiMgr stuff
#include "gui/GuiMgr.h"
#include "gui/Menu.h"
#include "gui/SubMenu.h"
#include "gui/SimpleCommand.h"

// declare the service factories for the GuiSvc
static SvcFactory<GuiSvc> a_factory;
const ISvcFactory& GuiSvcFactory = a_factory;


// ------------------------------------------------
// Implementation of the GuiSvc class
// ------------------------------------------------
/// Standard Constructor
GuiSvc::GuiSvc(const std::string& name,ISvcLocator* svc)
: Service(name,svc), m_guiMgr(0)
{
    
    // declare the properties and set defaults

    declareProperty ("size", m_size=-300);
    declareProperty( "pause_interval", m_pause_interval=0);
    declareProperty( "paused", m_paused=true);

    declareProperty("EvtMax", m_evtMax=0xFEEDBABE);


}



void GuiSvc::pause()
{
    m_guiMgr->pause();
}


// initialize
StatusCode GuiSvc::initialize () 
{
    using namespace gui;
    StatusCode  status =  Service::initialize ();

        
    // bind all of the properties for this service
    setProperties ();
    
    // open the message log
    MsgStream log( msgSvc(), name() );
 
    // create and start the gui manager
    m_guiMgr = new gui::GuiMgr(m_size, 2, m_pause_interval); 

    // use the incident service to register begin, end events
    IIncidentSvc* incsvc = 0;
    StatusCode sc = service ("IncidentSvc", incsvc, true);
    
    if( status.isFailure() ) return status;

    incsvc->addListener(this, "BeginEvent", 100);
    incsvc->addListener(this, "EndEvent", 0);

    // get the menu for additional buttons.
    gui::SubMenu& sub_menu= m_guiMgr->subMenu();

    sub_menu.addButton("set max event...",
        new SimpleCommand<GuiSvc>(this, &GuiSvc::queryEvtMax));

    m_guiMgr->menu().file_menu().addButton("set OutputLevel...",
        new SimpleCommand<GuiSvc>(this, &GuiSvc::queryOutputLevel));

    sub_menu.addButton("Quit Loop",
        new SimpleCommand<GuiSvc>(this, &GuiSvc::quit));

    sub_menu.addButton("Set Pause interval...", 
        new SimpleCommand<GuiMgr>(m_guiMgr, &GuiMgr::queryPause));

    m_guiMgr->menu().add(new MenuClient<GuiSvc>(this)); // schedule callback if exit button pressed



    sc = serviceLocator()->queryInterface(IID_IAppMgrUI, (void**)&m_appMgrUI);

    // get property from application manager
    if ( m_evtMax == (int)0xFEEDBABE )   {
        IProperty* propMgr=0;
        status = serviceLocator()->service("ApplicationMgr", propMgr );
        if( status.isFailure()) {
            log << MSG::ERROR << "Unable to locate PropertyManager Service" << endreq;
            return status;
        }
                
        IntegerProperty evtMax("EvtMax",0);
        status = propMgr->getProperty( &evtMax );
        if (status.isFailure()) return status;
        int max_event = evtMax.value();
        
        setProperty(evtMax);
    }
    
    //----------------------------------------------------------------
    // most of  the following cribbed from ToolSvc and ObjManager

    // look for a factory of an AlgTool that implements the IGuiTool interface:
    // if found, make one and call the special method 

    // Manager of the AlgTool Objects
    IObjManager* objManager=0;             

    // locate Object Manager to locate later the tools 
    status = serviceLocator()->service("ApplicationMgr", objManager );
    if( status.isFailure()) {
        log << MSG::ERROR << "Unable to locate ObjectManager Service" << endreq;
        return status;
    }

    
    IToolSvc* tsvc  =0;
    status = service( "ToolSvc", tsvc, true );
    if( status.isFailure() ) {
         log << MSG::ERROR << "Unable to locate Tool Service" << endreq;
        return status;
    }

    IToolFactory* toolfactory = 0;
    
    // search throught all objects (factories?)
    for(IObjManager::ObjIterator it = objManager->objBegin(); it !=objManager->objEnd(); ++ it){
        
        std::string tooltype= (*it)->ident();
        // is it a tool factory?
        const IFactory* factory = objManager->objFactory( tooltype );
        IFactory* fact = const_cast<IFactory*>(factory);
        status = fact->queryInterface( IID_IToolFactory, (void**)&toolfactory );
        if( status.isSuccess() ) {

            // yes: now see if the tool implements the IGuiTool interface
            IGuiTool* ireg;
            status = tsvc->retrieveTool(tooltype, ireg);
            if( status.isSuccess() ){
                log << MSG::DEBUG << "Initializing gui stuff in " << tooltype << endreq;
                ireg->initialize(m_guiMgr);
            }
            
        }

    }

    return StatusCode::SUCCESS;
}


void GuiSvc::queryEvtMax()
{

    m_guiMgr->menu().query("Enter new max event",& m_evtMax);
}

void GuiSvc::queryOutputLevel()
{
    IProperty* glastPropMgr=0;
    StatusCode status = service("MessageSvc", glastPropMgr, true);
    if( status.isFailure() ) return;
      
    IntegerProperty levelProp("OutputLevel",0);
    status = glastPropMgr->getProperty( &levelProp );
    if (status.isFailure()) return;

    int level = levelProp.value();
    m_guiMgr->menu().query("Enter new outputlevel(2=DEBUG, 3=INFO, 4=WARNING, 5=ERROR, 6=FATAL )",& level);

    levelProp =level;
    status = glastPropMgr->setProperty( levelProp );
    glastPropMgr->release();

}


void GuiSvc::quit() {
    m_evtMax=0;
    m_guiMgr->resume();
}

// handle "incidents"
void GuiSvc::handle(const Incident &inc)
{
#if 0 // don't need now that implementing Runable interface
    if( inc.type()=="BeginEvent")beginEvent();
    else if(inc.type()=="EndEvent")endEvent();
#endif
}


void GuiSvc::beginEvent() // should be called at the beginning of an event
{ 
    using namespace gui;
    static bool first=true;
    if( first) { first=false;

        // start of first event: startup the GUI
        m_guiMgr->setState( m_paused?  GuiMgr::PAUSED : GuiMgr::RUNNING);
        MsgStream log( msgSvc(), name() );

        log << MSG::INFO << "starting GuiSvc GUI in " 
            << (m_guiMgr->state()==GuiMgr::PAUSED? "PAUSED": "RUNNING" ) << " state" << endreq;
        m_guiMgr->menu().run();
    }
    m_guiMgr->begin_event();
}

void GuiSvc::endEvent()  // must be called at the end of an event to update, allow pause
{        
    m_guiMgr->end_event();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode GuiSvc::run(){
    if ( 0 != m_appMgrUI )    {
        // loop over the events
        int event= 0;
        while(event++ < m_evtMax) {
            beginEvent();
            if( m_appMgrUI->nextEvent(1).isFailure()) return StatusCode::FAILURE;
            endEvent();
        }
    }
    return StatusCode::SUCCESS;
    
}

// finalize
StatusCode GuiSvc::finalize ()
{
    StatusCode  status = StatusCode::SUCCESS;
       
    delete m_guiMgr;
    return status;
}
/// Query interface
StatusCode GuiSvc::queryInterface(const IID& riid, void** ppvInterface)  {
  if ( IID_IGuiSvc.versionMatch(riid) )  {
    *ppvInterface = (IGuiSvc*)this;
  }
  else if (IID_IRunable.versionMatch(riid) ) {
      *ppvInterface = (IRunable*)this;
  }
  else  {
    return Service::queryInterface(riid, ppvInterface);
  }
  addRef();
  return SUCCESS;
}

