//     $Id$
//  Author: Toby Burnett
//

#ifndef MENU_H
#define MENU_H

#include "gui/Command.h"
#include "gui/GUI.h"

#include <vector>
#include <map>

namespace gui {

class SubMenu;
class MenuBase;

class Menu {
// This class manages a hierarchy of <<Command>> instances, each associated with a 
// unique name. If there is an associated <<GUI>>, it will also enter each command into
// the GUI's menu bar
//
// When started, it establishes the File menu, and adds an Exit command to it when
// its run is called.
//
// The Menu also manages a table of keystroke commands
public:
    friend class SubMenu;

    Menu(GUI&);
    // constuctor with a reference to the GUI
    
    ~Menu();

    SubMenu& subMenu(const std::string& label);
    // Return a SubMenu (pull-down on a GUI menu bar)
      
    SubMenu& file_menu();
    // Return reference to the file menu sub menu object;

    void add_button(const std::string& label, Command* cmd);
    // Add top-level button (not kosher in Motif)

    void register_key(char key, Command* cmd);
    // Register a short-cut key with a  command

    bool strike(char c);
    // Execute the command, if any, associated with the key c. 
    // Return true if it has been registered, false otherwise

    void query(const std::string& ask, int* value);
    // Ask the user for an int   

    void query(const std::string& ask, double* value, int count=1);
    // Ask the user for a double (or a list)

    void query(const std::string& ask, float* value);
    // Ask the user for a float

    void query(const std::string& ask, std::string* value);
    // Ask the user for a string

    void run(bool paused=true);
    // register the "quit" command after any user file/IO comamnds
    // start the GUI, and optionally the message loop
    
    void check_messages();
    // check for message in the GUI (needed for long events if single threaded

    void add(MenuBase* client);
    // client uses to schedule a call to its finishSetup method
    
    static Menu* instance(Menu* m=0); 
	// deprecated method

    void quit();

private:

    class Node{
	//empty class used as tag
    };

    // following used by implementation
    void addCommand(const std::string& label, Command* command);
    GUI::Toggle* addToggle(const std::string& label, bool state, Command* set, Command* unset);
    Node* beginMenu(const std::string&, Node* subnode=0);
    void setMenu(Node* m);
    void endMenu();
    void addSeparator();

    static Menu* s_instance;

    GUI& m_gui; // reference to the GUI
    
    SubMenu*    m_fileMenu;
    // for adding file menu commands; set by the constructor

    typedef std::map<char, Command*, std::less<char> > KeyCommandMap;
    KeyCommandMap m_key_map;
    // list of keys and commands
    
    typedef    std::vector<MenuBase* > ClientList;
    ClientList	m_clients;
    // list of clients who want finishSetup to be called
    
    std::vector<SubMenu*> m_sub_menus;
    // list of sub menus to delete

    std::vector<Command*> m_commands;
    // list of commands to delete
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template< class T >
class SimpleDialogBox : public Command {
	// define a templated class that can be used as a command to query for a single
	// entity
public:
    SimpleDialogBox(Menu* m, std::string& q, T* t)
	: m_menu(m),
	m_question(q),
	m_t(t){}
    void execute(){m_menu->query(m_question, m_t);}
private:
    Menu* m_menu;
    std::string m_question;
    T* m_t;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class MenuBase {
// Virtual base class for the templated menu client class
public:
    MenuBase(){}
    virtual void finishSetup(){};
    virtual void quit(){};
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template< class T >
class MenuClient : public MenuBase {
    // Template class with the class of the client as the template parameter. 
    // If a client supports finishSetup  method, it can add
    // itself with the statement: "menu().add(new MenuClient<Client>(this));"
    // it will be called back just before starting if it implements finishSetup(),
    // and its quit() if the File/Exit button is pushed. 
public:
    MenuClient(T* t):m_t(t){}
	// Constructor taking the client as argument

    void finishSetup(){m_t->finishSetup();}
    void quit(){m_t->quit();}
	// Send the finishSetup() or quit() to the client
    friend class Menu;
private:
    T* m_t;
};

inline void Menu::add(MenuBase* client){ m_clients.push_back(client);}

} //namespace gui
#endif










