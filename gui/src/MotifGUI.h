// $Header$
//  Author: G. Barrand, T. Burnett
// Motif GUI definition.

#ifndef MotifGUI_H
#define MotifGUI_H

#include "gui/GUI.h"
#include "gui/SceneControl.h"

using namespace gui;

//class SceneControl;

class MotifGUI : public GUI {
// Derived GUI supporting the Motif interface
public:
    MotifGUI(const char* app="Motif", const char* title=0);

    void addToMenu(const char* title, Command* command);
    virtual GUI::Menu* beginPullDownMenu(const char* title, Menu* m=0);

    void restorePullDownMenu(GUI::Menu* m);
    GUI::Toggle* addToggleToMenu(const char* title, bool state, Command* set, Command* unset);
    void endPullDownMenu();
    virtual void menuSeparator();
    // hand mmenus

	gui::SceneControl* graphicsWindow(float size, int inital_view );
	std::ostream* textWindow(const char* name);
    // create, return handle to associated window

    void setTitle(const char* newtitle);

    void processMessages();
    // call to process waiting messages

    void run(int pause_interval);
    // start the window system

    void quit();

    char * askForFileName(const char* a, const char* b, const char* c);
    char* askUser(const char* promptString, const char* defaultString);
    void inform(const char* msg);

    const char * nameOf()const{return "MotifGUi";}


private:

    SceneControl* m_scene;
    class Toggle;

};
#endif

