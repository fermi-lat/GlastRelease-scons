//     $Id$
//  Author: Toby Burnett
//

#ifndef SUBMENU_H
#define SUBMENU_H
#include <string>
#include "gui/Command.h"
#include "gui/Menu.h"

namespace gui{

class SubMenu {
    // Handles commands in a nested pull-down menu, either from the menu bar, or of a submenu
public:
    void addButton(const std::string& label, Command* command);
    // add a simple button, with name label. Calls the command when struck

    GUI::Toggle* addToggle(const std::string& label, bool state, Command* set, Command* unset);
    // add a toggle button, initialy set if state is true. the comands set and unset are
    // called depending on the current state 
    // return pointer to object that can be used to access the state

    void addSeparator();
    // add a separator

    SubMenu& subMenu(const std::string& label);
    // Return a nested SubMenu for this submenu

    const Menu::Node* tag()const{ return m_tag;} 
    // return tag assigned by the windowing system, a pointer to something that the 
    // client must cast

private:
    friend class Menu; // only Menu has access to the constructor
    SubMenu(Menu& m,    const std::string& label);
    SubMenu(Menu& m, SubMenu* s, const std::string& label);

    Menu& m_menu;       // top level menu bar
    SubMenu* m_subMenu; // parent submenu, null if parent is menu
    Menu::Node* m_tag;  // used to communicate with real gui
};
} //gui

#endif //SUBMENU
