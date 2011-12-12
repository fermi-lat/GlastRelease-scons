/** 
* @file RootMenu.cxx
* @brief declaration and implementation of class RootMenu
*
*  $Header$
*/

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"

#include "GuiSvc/IGuiTool.h"
#include "RootIo/IRootIoSvc.h"

#include "gui/GuiMgr.h"
#include "gui/SubMenu.h"
#include "gui/GuiMgr.h"
#include "gui/Command.h"
#include "gui/SimpleCommand.h"

/** 
* \class RootMenu
*
* \brief Setup menu for requesting a random event from Root reading
* \author Heather Kelly
* 
*  

* $Header$
*/

class RootMenu : public AlgTool, virtual public IGuiTool {
 public:
     
   RootMenu( const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~RootMenu() { }
    
    /// implement to define gui elements: will be called from GuiSvc
    StatusCode initialize(gui::GuiMgr*);


    // call back for created button
    void setEventIndex();
    void setRunEventPair();

    /// call back for the Menu.
    void finishSetup(); 
    void quit(){};

    int index()const {return m_index;}
    std::pair<int, int> runEvent()const {return m_runEventPair; }
    IRootIoSvc* rootIoService()const{return m_rootioSvc;}

private:

    // Data Members
    gui::GuiMgr* m_guiMgr;
    IRootIoSvc* m_rootioSvc;
	int m_index;
	std::pair<int,int> m_runEventPair;

};


// Static factory for instantiation of algtool objects
//static ToolFactory<RootMenu> s_factory;
//const IToolFactory& RootMenuFactory = s_factory;
DECLARE_TOOL_FACTORY(RootMenu);

// Standard Constructor
RootMenu::RootMenu(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
     : AlgTool( type, name, parent ) 
    , m_index(0)
{
    
	m_runEventPair = std::pair<int, int>(0,0);
    // Declare additional interface
    declareInterface<IGuiTool>(this);
        
}


StatusCode RootMenu::initialize(gui::GuiMgr* guiMgr) 
{
    m_guiMgr = guiMgr;

    // schedule callback to add the menu, so sure that Flux is initialized
    m_guiMgr->menu().add(new gui::MenuClient<RootMenu>(this)); 

    return StatusCode::SUCCESS;
} 


void RootMenu::finishSetup()
{
    using namespace gui;

    MsgStream log( msgSvc(), name() );

    IService* theSvc;
    StatusCode sc = serviceLocator()->getService( "RootIoSvc", theSvc, true );
    if ( sc.isSuccess() ) {
      sc = theSvc->queryInterface(IRootIoSvc::interfaceID(), (void**)&m_rootioSvc);
    }else {
        log << MSG::ERROR << "Could not get RootIoSvc" << endreq;
        return;
    }

    // now get the filemenu and add a source button to it.
    gui::SubMenu& filemenu = m_guiMgr->menu().file_menu();

    filemenu.addButton("set Root Index...", 
		new SimpleCommand<RootMenu>(this, &RootMenu::setEventIndex));

	// save this for later...once BuildIndex trouble is sorted out
    filemenu.addButton("set Root Run and Event...", 
		new SimpleCommand<RootMenu>(this, &RootMenu::setRunEventPair));

}

void RootMenu::setEventIndex()
{
    // purpose: callback for buttons to set a root index
    // Loop until we find a good index
    bool DONE = false;
    while (!DONE) {
        m_index = -99;
        m_guiMgr->menu().query("Enter Root Index", &m_index);
        if (m_index == -99) return; 
        DONE = m_rootioSvc->setIndex(m_index);
        
    }
}

void RootMenu::setRunEventPair() {
    int run, evt;
    bool DONE = false;
    // Loop until we find a good run/event id pair
    while (!DONE) {
        run = -99;
        evt = -99;
        m_guiMgr->menu().query("Enter Root Run Number", &run);
        m_guiMgr->menu().query("Enter Root Event Number", &evt);
        if ((run == -99) || (evt == -99)) return;
        m_runEventPair = std::pair<int, int>(run, evt);
        // Make sure Index is set to -1
        m_rootioSvc->setIndex(-1);
        DONE = m_rootioSvc->setRunEventPair(m_runEventPair);
    }
}
