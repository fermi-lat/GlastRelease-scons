//     $Id$
//  Author: Toby Burnett
//

#ifndef GUIMGR_H
#define GUIMGR_H

#include "gui/SubMenu.h"
#include "gui/PrintControl.h"
#include "gui/DisplayControl.h"
#include "gui/GUI.h"

#include <vector>
#include <iostream>

namespace gui {
        
    /** Sets up a visualization and control environment.
     GuiMgr also has a command source, encapsulated by a Menu and the
     output devices: DisplayControl and PrintControl.
    */
class GuiMgr {
public:
    
    /** Constructor: size will be the initial size of the display viewport
       which should enclose subsequent volumes
     initial view is the starting view number
     */
    GuiMgr( float size=2000, int initial_view=1, int pause_interval=-1 );
    
    ~GuiMgr();
    
    //! should be called at start of event
    void                begin_event();

    //! must be called at the end of an event to allow regaining control 
    void                end_event();

    //! if called from a user program, will update and refresh the display,
    //! then, if GuiMgr is in a paused state, go into a GUI wait loop until CR or spacebar. 
    //! The optional title bar text will be put on the title.  Call with empty string to
    //! keep same display.
    void                break_point(const std::string& title_bar_text = "[user pause]");

    //! return true if end of run, false to continue
    bool                done() const;
 
    
    //! Access to singleton instance
    static GuiMgr*        instance();
    
    //! Returns the GUI
    GUI& 		gui()const;
    
    //! Returns the DisplayControl
    DisplayControl& 	display()const;
    
    //! Returns the PrintControl
    PrintControl& 	printer()const;
    
    //! Returns the Menu
    Menu&	   	menu()const;

    //! returns the eventloop submenu
    SubMenu&            subMenu()const;
    
    //! access to world size
    float		size()const;
                
    // these needed to implement MenuClient interface
    void quit();
    void finishSetup(){};  // dummy


    // these connected to keys, buttons
    void pause();  // set status to pause on next loop (or single step)
    void resume(); // set status to skip next loop (or resume if paused)
    void stop_loop();
    void queryPause();
    enum State {INITIAL, RUNNING, PAUSED, DONE};
    State state()const;
    void setState(State s);

private: 

    GUI*    		m_gui;
    Menu*		m_menu;
    DisplayControl* 	m_display;
    PrintControl*   	m_print;
    float		m_size;
    SubMenu*            m_sub_menu;
    int			m_pause_interval;
    State               m_state;
    void setup(float size, int initial_view, int pause_interval);
    static GuiMgr* 	s_instance;
};

// inline trivial access functions
inline float            GuiMgr::size()const   {return m_size;}
inline GUI& 		GuiMgr::gui()const    {return *m_gui; }
inline DisplayControl&  GuiMgr::display()const{return *m_display;}
inline PrintControl& 	GuiMgr::printer()const{return *m_print;}
inline Menu& 		GuiMgr::menu()const   {return *m_menu;}
inline SubMenu&         GuiMgr::subMenu()const{return *m_sub_menu;}
inline bool             GuiMgr::done()const   {return m_state==DONE;}
inline gui::GuiMgr::State    GuiMgr::state()const  {return m_state;} 
inline void             GuiMgr::setState(gui::GuiMgr::State state) {m_state=state;}
 
// inline the constructor to force the GUI::createGUI call into the user's 
// application

inline GuiMgr::GuiMgr( float size, int initial_view, int pause_interval)
    :m_gui(GUI::createGUI("GuiMgr", "GuiMgr"))
{ 
    setup(size,initial_view,pause_interval);
}
} // namespace gui


#endif

