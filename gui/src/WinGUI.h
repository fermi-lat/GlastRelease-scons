// $Header$
//  Author: T. Burnett
//
// Win GUI definition.

#ifndef WinGUI_H
#define WinGUI_H

#if !defined (WIN32) && !defined(_WINDOWS)
#  error "This code only compiles for the Windows NT or Window95"
#endif

#define STRICT
#include <windows.h>
#undef min
#undef max
#include <vector>

#include "gui/GUI.h"
#include "gui/SceneControl.h"

using namespace gui;

#pragma warning(disable:4786)


class WinScene;
class WinGUIostream;

class WinGUI : public GUI {
    // Windows 95/NT implementation of GUI functions
public:
    WinGUI(const char* application="Win32", const char* title=0);
    ~WinGUI();
    
    void addToMenu(const char* title, Command* command);
    Menu* beginPullDownMenu(const char*, Menu* );
    GUI::Toggle* addToggleToMenu(const char* title, bool state, Command* set, Command* unset);
    void restorePullDownMenu(Menu*);
    void endPullDownMenu();
    void menuSeparator();
    // handle menus
    
    gui::SceneControl* graphicsWindow(float size, int initial_view);
    std::ostream* textWindow(const char* name);
    // create, return handle to associated window
    
    void start();
    // start the window system
        
    void run(int pause_interval);
    // start the message loop; return only when exit message

    void setTitle(const char* newtitle);
    void quit();
    
    void processMessages();
    
    char * askForFileName(const char* a, const char* b, const char* c);
    char* askUser(const char* promptString, const char* defaultString);
    void inform(const char* msg);
    
    const char* nameOf()const{return "WinGUI";}
    
    void execute(WPARAM);
    // execute the given command, needed by the Windows message handling procedure
    
    // following are static and public for communication with windows routines
    static HINSTANCE s_instance; // handle to the application instance. Set by windows main 
    static WinGUI* s_gui;  //pointer to instance 
    static HWND s_hwnd;     // windows handle to the main window
    static WinScene* s_graphics_window; // primary (3d display) window
    static WinScene* s_2d_window; // secondary (for plots) window
    static bool s_quitting;
    static WinGUIostream* s_text_window;

private:
    enum{ cmdOffset=100}; // constant for command indeces

    HMENU m_currentMenu;
    HMENU m_menuBar;

    UINT m_menuId;
    
    typedef std::vector<Command* > CommandList;
    CommandList m_commands;
    
    std::ostream* m_text_window;
    SceneControl* m_graphics_window;

    CommandList m_owned_commands; // separate list for deleting the extra toggle commands
    
    class Toggle; // used in implementation
};
#endif

