//     $Id$
//  Author: Toby Burnett
//
//
// Implement the PrintControl class 
//
// When this class has been set up, new print commands can be added to the menu
// by:      
//      PrintControl::instance()->addItem("my printout", myCommand);
//
// where
//    myCommand is an instance of a subclass of Command that will output to
// the std::ostream* PrintControl::instance()->out(); 


#include "gui/SubMenu.h"
#include "gui/PrintControl.h"
#include "gui/SimpleCommand.h"


#include <cassert>

namespace gui {
//////////////////////////////////////////////////////////////////////////////
//   Constructor: set up menu, get the GUI's text window to print to
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PrintControl::PrintControl(Menu& menu, std::ostream*out )
    : m_print(out), 
      m_menu(menu)
{
    setup();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PrintControl::~PrintControl()
{
    PrintList::iterator it = m_print_commands.begin(); 
    while (it !=  m_print_commands.end()) 
	delete (*it++);

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PrintCommand::~PrintCommand()
{ 
    delete m_cmd; 
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PrintControl::setup()
{
    s_instance = this;
    disable();
    // create Print pull-down menu
    m_sub_menu = &m_menu.subMenu("Print");
    // schedule call back to finishSetup method
    m_menu.add(new MenuClient<PrintControl>(this));

}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// file scope class that implements the command, followed by a flush

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PrintControl::add(Command* cmd, const std::string& title,  bool enabled)
{

    PrintCommand* prt = new PrintCommand(cmd, this, enabled);
    m_print_commands.push_back(prt);

    m_sub_menu->addToggle(title, enabled  
	,new SimpleCommand<PrintCommand>(prt, &PrintCommand::enable)
	,new SimpleCommand<PrintCommand>(prt, &PrintCommand::disable) );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PrintControl::enable(){    m_enabled = true;}
void PrintControl::disable(){    m_enabled = false;}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void 
PrintControl::printAll()
{
    if(m_enabled) printNow();
}
void
PrintControl::printNow()
{
    for( PrintList::iterator it = m_print_commands.begin(); it !=  m_print_commands.end(); ++it) 
	(*it)->execute();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void
PrintControl::finishSetup()
{

    m_sub_menu->addSeparator();
    m_sub_menu->addToggle("print after event", true, 
	new SimpleCommand<PrintControl>(this,&PrintControl::enable),
	new SimpleCommand<PrintControl>(this,&PrintControl::disable));
    m_sub_menu->addButton("print now", new SimpleCommand<PrintControl>(this,&PrintControl::printNow) );

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PrintControl::clear(){m_print->clear();}

void PrintControl::flush(){m_print->flush();}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//   stuff for singleton 
PrintControl* PrintControl::s_instance = 0;

PrintControl* PrintControl::instance()
{
    assert( s_instance!=0)	;
    return s_instance;
}

} // namespace gui
