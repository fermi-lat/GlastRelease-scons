//     $Id$
//  Author: Toby Burnett
//
// Define the singleton PrintControl class 
//
#ifndef PRINTCONTROL_H
#define PRINTCONTROL_H

#include "gui/Menu.h"

#include <vector>
#include <iostream>

namespace gui {

class PrintCommand; // defer to define in .cxx since private

class PrintControl  {
    // When this class has been set up, new print commands can be added to the menu
    // by:
    //      
    //      PrintControl::instance()->addPrinter("my printout", myCommand);
    //
    // where myCommand is an instance of a subclass of Command that
    // will output to the std::ostream* PrintControl::instance()->out();
public:
    PrintControl(Menu& menu, std::ostream* out);
    // constructor--specify a stream to print to

    virtual ~PrintControl(); 
    // dtor 

    void addPrinter(const std::string& label, Command * cmd, bool enabled=0);
    void add(Command * cmd, const std::string& label,  bool enabled=0);
    // adds a menu item, associated with any command, presumably a printer

    void printAll();
    // print if enabled -- intended for after event

    void enable();

    void disable();
    
    void printNow();
    // print anyway

    void clear();
    // clear the display device (may not be implemented)

    void flush();
    // pass on to implementor

    std::ostream* out(){return m_print;}
    // get at the ostream

    void setOut(std::ostream* o){m_print = o;}
    // replace the ostream object

    // these implement the MenuClient interface
    void finishSetup();    // called by the Menu object after user menu commands executed
    void quit(){};  // dummy for now

    const char* nameOf()const {return "PrintControl";}
    
private:
    void setup();
    
    std::ostream* m_print;
    // the stream to print to
        
    Menu& m_menu;
    // reference to our menu
    
    SubMenu* m_sub_menu;
    // pointer to print sub menu

    bool m_enabled;
    // determines if printing occurs after event

    typedef std::vector< PrintCommand* > PrintList;

    PrintList m_print_commands;


};

/** Template class that makes it easy to add a printer. Just do,

   PrintControl::instance()->addPrinter("My stuff",new Printer_T<Me>(this));

 and implement void printOn(ostream& out)const;
 */
template<class T>
class Printer_T : public Command {
public:
    Printer_T(const T* t, std::ostream* out=&std::cout):m_t(t),m_out(out){}
    void execute(){m_t->printOn(*m_out);}
private:
    const T* m_t;
    std::ostream* m_out;
};


class PrintCommand : public Command {
    // private class for PrintControl--needs to be here
    friend class PrintControl;
    PrintCommand( Command* cmd, PrintControl* printer, bool enabled )
	:m_cmd(cmd), m_printer(printer), m_enabled(enabled){}
    ~PrintCommand() ;
    void execute() {
	if(m_enabled){ 
	    m_cmd->execute();	
	    m_printer->flush();
	}
    }
    void enable(){m_enabled=true;}
    void disable(){m_enabled=false;}
    Command* m_cmd;
    PrintControl* m_printer;
    bool m_enabled;
};

// here for temporary 
inline void PrintControl::addPrinter(const std::string& label, Command * cmd, bool enabled)
{
    add(cmd, label, enabled);
}

} // namespace gui
#endif

