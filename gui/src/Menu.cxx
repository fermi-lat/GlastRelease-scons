//     $Id$
//  Author: Toby Burnett
//
//   Menu class implementations, now just pass-through to GUI


//#include "gui/Menu.h"
#include "gui/SubMenu.h"
#include "gui/SimpleCommand.h"
#include <strstream>
#include <string>
#include <cstdio> // for sprintf

namespace gui {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Menu::Menu(GUI& gui)
    :m_gui(gui)
{
    // create the pull-down File menu, save a pointer for adding commands later
    m_fileMenu = &subMenu("File");
    endMenu();
    s_instance = this;
}
Menu::~Menu()
{
    std::vector<SubMenu*>::iterator psm = m_sub_menus.begin();
    while (psm != m_sub_menus.end())
	delete (*psm++);

    ClientList::iterator it = m_clients.begin();
    while (it != m_clients.end())
	delete (*it++);
    
    std::vector<Command*>::iterator pcmd = m_commands.begin();

    while( pcmd != m_commands.end()) {
        delete (*pcmd++);
    }
	
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Menu* Menu::s_instance = 0;

Menu* Menu::instance(Menu* m)
{
    if(m) s_instance = m;
    return s_instance;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Menu::run(bool paused)
// start up the application
{
    // have all clients finish their setup
    ClientList::iterator it = m_clients.begin();
    while (it != m_clients.end())
	(*it++)->finishSetup();
    
    // add quit to the file menu
    m_fileMenu->addSeparator();
    m_fileMenu->addButton("Exit",   new SimpleCommand<Menu>(this, &Menu::quit));

    // run the gui message loop, or whatever it does
    if( !paused ) m_gui.start(); // if not paused
    else m_gui.run();
}
void Menu::quit()
{
    // have all clients quit
    ClientList::iterator it = m_clients.begin();
    while (it != m_clients.end()) (*it++)->quit();

    // then shut down message loop
    m_gui.quit();
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Menu::check_messages()
{
    m_gui.processMessages();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Menu::add_button(const std::string& title, Command* command)
{
    m_gui.addToMenu( title.c_str(),  command);
    m_commands.push_back(command);  // save to maybe delete

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Menu::addCommand(const std::string& title, Command* command)
{
    m_gui.addToMenu( title.c_str(),  command);

    m_commands.push_back(command);  // save to maybe delete
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GUI::Toggle* Menu::addToggle(const std::string& title, bool state, Command* set, Command* unset)
{
    GUI::Toggle* t =m_gui.addToggleToMenu(title.c_str(), state, set, unset);
    if( state ) set->execute(); else unset->execute();
    m_commands.push_back(set);
    m_commands.push_back(unset);
    return t;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Menu::Node* Menu::beginMenu(const std::string& name, Menu::Node* subnode)
{
    return (Menu::Node*)m_gui.beginPullDownMenu(name.c_str(), (GUI::Menu*) subnode);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Menu::setMenu(Node* m)
{
    m_gui.restorePullDownMenu( (GUI::Menu*) m);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SubMenu& Menu::file_menu()
{
    return *m_fileMenu;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Menu::addSeparator()
{
    m_gui.menuSeparator();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Menu::endMenu()
{
    m_gui.endPullDownMenu();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  shortcut key implementations
void Menu::register_key(char key, Command* cmd)
{
    m_key_map[key]=cmd;
}
bool Menu::strike(char key)
{
    KeyCommandMap::const_iterator it = m_key_map.find(key);
    if( it != m_key_map.end() ) {
        (*it).second->execute();
	return true;
    }
    return false;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Menu::query(const std::string& ask, int* value)
{
    char temp[20]; sprintf(temp,"%i", *value);
    char* answer = m_gui.askUser(ask.c_str(),temp ) ;
    if( answer && answer[0] ) 
	    *value =atoi(answer) ; 
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Menu::query(const std::string& ask, double* value, int n)
{
    std::ostrstream temp;
    for(int i=0; i<n; ++i) {
        temp << value[i] ;
        if( i<n-1) temp << ", ";
    }
    temp << '\0';
    std::string answer(m_gui.askUser(ask.c_str(), temp.str() ) );
#ifdef _MSC_VER // not defined otherwise?
    temp.rdbuf()->freeze(false); // allow delete
#endif
    if( answer.size()==0 ) return; // check for cancel
    answer += ","; // for convenience 
    //tokenize
    int p = 0; 
    for(int j=0; j< n; ++j) {
        int q = answer.substr(p).find_first_of(",");
        if( q <0 ) break;
        if( q >0 ) value[j] = atof(answer.substr(p,q).c_str());
        p += q+1; // skip to next field
    }


}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Menu::query(const std::string& ask, float* value)
{
    char temp[20]; sprintf(temp,"%f", *value);
    char* answer = m_gui.askUser(ask.c_str(),temp ) ;
    if( answer && answer[0] ) 
	    *value =float(atof(answer)) ; 
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Menu::query(const std::string& ask, std::string* value)
{
    std::string answer = m_gui.askUser(ask.c_str(),value->c_str() ) ;
    if( answer.length() >0 ) 
	    *value =answer ; 
}

} // namespace gui
