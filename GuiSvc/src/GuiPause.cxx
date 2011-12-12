/** 
* @file GuiPause.cxx
* @brief force a pause if the GUI is running
*/
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Property.h"

// if use the gui
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"
#include "gui/DisplayControl.h"
#include "gui/SimpleCommand.h"

#include <string>
#include <iostream>


/** @class GuiPause
* @brief Algorithm to pause execution to  user to examine an event more closely. 
*
* $Header$
<br> Example jobOptions fragment:
<pre>
</pre>
*/
class GuiPause : public Algorithm {
public:
    GuiPause(const std::string& name, ISvcLocator* pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

    void enable(){
            m_enabled=true;
            m_guiSvc->guiMgr()->setState(gui::GuiMgr::SKIPPING);
    }
    void disable(){
            m_enabled=false;
            m_guiSvc->guiMgr()->setState(gui::GuiMgr::PAUSED);
    }
private: 

    //! access to the Gui Service for display of 3-d objects
    IGuiSvc*    m_guiSvc;

    bool m_enabled;
};

// necessary to define a Factory for this algorithm
// expect that the xxx_load.cxx file contains a call     
//     DLL_DECL_ALGORITHM( GuiPause );

//static const AlgFactory<GuiPause>  Factory;
//const IAlgFactory& GuiPauseFactory = Factory;
DECLARE_ALGORITHM_FACTORY(GuiPause);

GuiPause::GuiPause(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator), m_guiSvc(0) , m_enabled(false)
{
}

StatusCode GuiPause::initialize()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    // get the Gui service
    if (service("GuiSvc", m_guiSvc).isFailure ()){
        log << MSG::WARNING << "GuiSvc not found" << endreq;
    }else{
        // enable skiping
        gui::DisplayControl& display = m_guiSvc->guiMgr()->display();
        display.menu().addToggle("enable skip", false,
            new gui::SimpleCommand<GuiPause>(this,&GuiPause::enable),
            new gui::SimpleCommand<GuiPause>(this,&GuiPause::disable));
    }
    return sc;
}

StatusCode GuiPause::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    if( m_enabled){
    // have it stop on on this event
         if(m_guiSvc!=0) m_guiSvc->guiMgr()->setState(gui::GuiMgr::INTERRUPT);
    }
    return sc;
}

StatusCode GuiPause::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    return sc;
}

