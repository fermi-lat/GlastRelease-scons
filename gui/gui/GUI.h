//     $Id$
// Project: Atlas top-level prototype
//  Author: Toby Burnett

#ifndef __GUI_H
#define __GUI_H

#include "gui/Command.h"

namespace gui { class SceneControl; }

#include <iostream>

namespace gui {

class GUI 
// Abstract base class which defines a factory method, createGUI for
// the "real" GUI.
// Defines methods to set the window title, to define a menu bar with pull-down command menus,
// (call-backs to Command objects), to use a graphical display (manipulated by an appropriate
// SceneControl object), and to direct printout to a separate window.
{
public:
    virtual ~GUI();

    virtual void addToMenu(const char* title, Command* command);
    // Add a button to the current menu, either the menu bar, or a pull-down menu

    
    class Toggle : public Command 
    // Special abstract toggle command that must be implemented by the specific GUI subclass
    // set and unset methods allow control of state outside of the windowing system, so that
    // a button can control other buttons, setting or unsetting them independent of their 
    // current state.
    {
    public:
        void set()  { if(!state()) execute();}
        void unset(){ if( state()) execute();}
        virtual bool state()const=0;
    };
    
    virtual Toggle* addToggleToMenu(const char* title, bool state, Command* set, Command* unset);
    // Add a toggle button to the current menu. state is the initial state (on or off)
    // return pointer to command that is created to respond to mouse click

    class Menu; // empty class to tag menu bar or sub menu

    virtual Menu* beginPullDownMenu(const char* title, Menu* m=0);
    // start adding commands to a pull-down menu, which may be a submenu of a previous one.
    // If zero, add to menubar.
    // Return tag

    virtual void restorePullDownMenu(Menu* m);
    // Restore a previously-defined pull-down menu

    virtual void endPullDownMenu();
    // Stop adding commands to a pull-down menu and revert to the menu bar

    virtual void menuSeparator();
    // Insert a separator bar in the current pull-down menu

    virtual gui::SceneControl* graphicsWindow(float size=100, int initial_view=1)=0;
    // Create and return a SceneControl object that can create gui graphics in the graphical
    // portion of the window.

    virtual std::ostream* textWindow(const char*)=0;
    // Create and return an ostream object. (Which may just be cout for unix.)

    virtual void start();
    virtual void run(int pause_interval=-1);
    // Start the window system. For a real GUI, this is a message loop.

    virtual void quit();
    // End the window system. Typically bound to a "Exit" button.

    virtual void setTitle(const char* newtitle);
    // Set the window title

    virtual void processMessages();
    // Process waiting messages (needed if single thread)

    virtual char * askForFileName(const char* a, const char* b, const char* c);
    // File name dialog box.

    virtual char * askUser(const char* promptString,const char* defaultString);
    // Simple dialog box that sets a string

    virtual void inform(const char* msg);
    // Dialog box with a message.

    static GUI* createGUI(const char* appname=0, const char* title=0);
    // Return a GUI. Expect subclass to implement. Prevents more that one GUI implementation
    // in an executable

    static GUI* instance(GUI* i=0);
    static bool running;  // use to check that GUI has been started

protected:
    GUI(const char* title);  // abstract class
    const char* _winTitle;   // window title
private:
    static GUI* s_instance;  // public to allow simple setup

};

inline GUI::GUI(const char* title):_winTitle(title){}  // abstract class
inline GUI::~GUI(){}
inline void    GUI::addToMenu(const char* , Command* ){}
inline GUI::Toggle*   GUI::addToggleToMenu(const char*, bool, Command*, Command* ){return 0;}
inline GUI::Menu*   GUI::beginPullDownMenu(const char*, GUI::Menu* m){return 0;}
inline void    GUI::restorePullDownMenu(GUI::Menu* ){}
inline void    GUI::endPullDownMenu(){}
inline void    GUI::menuSeparator(){}
inline void    GUI::start(){}
inline void    GUI::run(int){}
inline void    GUI::quit(){}
inline void    GUI::setTitle(const char*){}
inline void    GUI::processMessages(){}
inline char *  GUI::askForFileName(const char*, const char*, const char* ){return 0;}
inline char *  GUI::askUser(const char* ,const char* ){return 0;}
inline void    GUI::inform(const char* ){}


} // namespace gui

#endif

