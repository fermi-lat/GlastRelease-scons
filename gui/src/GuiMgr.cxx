//     $Id$
//  Author: Toby Burnett
//
//  Implement GuiMgr

#include "gui/GuiMgr.h"
#include "gui/SimpleCommand.h"

#ifdef WIN32
#include <sstream>
#else // gcc < 2.9.5.3 dosen't know <sstream>,stringstream
#include <strstream>
#endif

#include <iostream>


namespace gui {


GuiMgr* GuiMgr::s_instance = 0;
GuiMgr* GuiMgr::instance(){ if (s_instance) return s_instance; else return new GuiMgr;} 

static std::string title; // keep title for refenece 

void GuiMgr::setup(float size, int inital_view, int pause_interval )
{
    s_instance = this;
    
    m_state = PAUSED;
    
    m_menu = new Menu(*m_gui);
    m_display =new DisplayControl(menu(),  m_gui->graphicsWindow(0.5*size, inital_view) ) ;
    m_print   =new PrintControl(menu(), m_gui->textWindow("output") );
    
    m_size = size ;
    
    // Define submenu on menu bar (which is available for more buttons) 
    m_sub_menu =&menu().subMenu("EventLoop");
    
    Command* resume=new SimpleCommand<GuiMgr>(this,&GuiMgr::resume);
    m_sub_menu->addButton("Start/Resume (cr)", resume);
    
    Command* pause = new SimpleCommand<GuiMgr>(this,&GuiMgr::pause);
    m_sub_menu->addButton("pause/step (space)",pause);
    
    m_pause_interval=pause_interval; 
    m_sub_menu->addButton("set pause interval...", new SimpleCommand<GuiMgr>(this,&GuiMgr::queryPause));
    
    // now set up some useful key-stroke shortcuts
    menu().register_key(' ', pause);    // space bar -- single event (or pause run)
    menu().register_key('\r',resume);   // cr -- start or resume
    
    menu().add(new MenuClient<GuiMgr>(this)); // schedule callback if exit button pressed
}

void GuiMgr::queryPause()
{
    menu().query("Enter new pause interval in ms (0: no pause, -1 infinite)", &m_pause_interval);	
}



GuiMgr::~GuiMgr()
{    
    delete  m_gui;
    delete  m_display;
    delete  m_print;
    delete  m_menu;
}



void GuiMgr::quit() // called back when File/Exit button pushed
{
    stop_loop();
    m_state = DONE;   
}


void GuiMgr::begin_event() // should be called at the beginning of an event
{ 
    static int n=0;
    display().clear();
#ifdef WIN32
    std::stringstream label; label << " Event " << (++n) << '\0';
#else // gcc < 2.9.5.3 dosen't know <sstream>,stringstream
    std::strstream label; label << " Event " << (++n) << '\0';
#endif
    title = std::string(label.str());
    display().setTitle(title); 
}

void GuiMgr::end_event()  // must be called at the end of an event to update, allow pause
{        
    display().update();
    printer().printAll();
    gui().processMessages();
    
    if( ! done() ) {
        // start a new message loop: process messages until GUI::running set false
        gui().run(m_state == RUNNING ? m_pause_interval : -1);
    }
}

void GuiMgr::break_point(const std::string& title_bar_text)
{
    if( !title_bar_text.empty() ) display().setTitle(title+title_bar_text);
    end_event();
    if( !title_bar_text.empty() ){ display().setTitle(title); gui().processMessages();}
}

void GuiMgr::pause() // called by GUI message (space bar or button)
{   
    switch  (m_state) {
    case RUNNING: m_state=PAUSED;  stop_loop(); break;
    case PAUSED:  stop_loop(); break;
    }
}

void GuiMgr::resume() // called by GUI message (CR key or button )
{
    switch  (m_state) {
    case RUNNING: m_state=PAUSED;  stop_loop(); break;
    case PAUSED:  m_state=RUNNING; stop_loop(); break;
    }
}
void GuiMgr::stop_loop() // stop the loop in the end_event method
{
    GUI::running = false;
}



} // namespace gui








