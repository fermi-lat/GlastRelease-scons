// $Header$

#include "gui/SubMenu.h"

#include <string>
using namespace gui;

//==============================================================================
//        SubMenu implementation              
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SubMenu& Menu::subMenu(const std::string& name)
{
    SubMenu* sm = new SubMenu(*this, name);
    m_sub_menus.push_back(sm);
    return * sm;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SubMenu::SubMenu(Menu& m, const std::string& name)
    :m_menu(m), m_subMenu(0),
     m_tag(m.beginMenu(name))
{
	m.endMenu();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SubMenu::SubMenu(Menu& m, SubMenu* s, const std::string& name)
    :m_menu(m), 
    m_subMenu(s),
     m_tag(m.beginMenu( name, s->m_tag))
{
	m.endMenu();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SubMenu::addButton(const std::string& name, Command* cmd)
{
    m_menu.setMenu(m_tag);
    m_menu.addCommand(name,cmd);
    m_menu.endMenu();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SubMenu::addSeparator()
{
    m_menu.setMenu(m_tag);
    m_menu.addSeparator();
    m_menu.endMenu();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GUI::Toggle* SubMenu::addToggle(const std::string& title, bool state, Command* set, Command* unset)
{
    m_menu.setMenu(m_tag);
    GUI::Toggle* t = m_menu.addToggle(title,state,set,unset);
    m_menu.endMenu();
    return t;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SubMenu& SubMenu::subMenu(const std::string& label)
{
    // creating a submenu of this menu, but looks the same to the menu
    
    SubMenu* sm = new SubMenu(m_menu, this, label);
    m_menu.m_sub_menus.push_back(sm);
    return * sm;
}
