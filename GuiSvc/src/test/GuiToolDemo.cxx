#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"

#include "GuiSvc/IGuiTool.h"

#include "gui/GuiMgr.h"
#include "gui/SubMenu.h"
#include "gui/GuiMgr.h"
#include "gui/Command.h"

class TestMenu : public AlgTool, virtual public IGuiTool {
 public:
     
   TestMenu( const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TestMenu() { }
    
    /// implement to define gui elements: will be called from GuiSvc
    StatusCode initialize(gui::GuiMgr*);


    // call back for created button
    void setSource(std::string newsource);

    /// call back for the Menu.
    void finishSetup(); 
    void quit(){};


private:

    // Data Members
    gui::GuiMgr* m_guiMgr;

};


// Static factory for instantiation of algtool objects
static ToolFactory<TestMenu> s_factory;
const IToolFactory& TestMenuFactory = s_factory;

// Standard Constructor
TestMenu::TestMenu(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : AlgTool( type, name, parent ) {
    
    // Declare additional interface
    declareInterface<IGuiTool>(this);
        
}


StatusCode TestMenu::initialize(gui::GuiMgr* guiMgr) 
{
    m_guiMgr = guiMgr;

    // schedule callback to add the menu
    m_guiMgr->menu().add(new gui::MenuClient<TestMenu>(this)); 


    return StatusCode::SUCCESS;
} 


void TestMenu::finishSetup()
{

    MsgStream log( msgSvc(), name() );

    
    //   ----------------------------------
    class TestCommand : public gui::Command { 
    public:
        TestCommand(gui::GUI& gui):m_gui(gui){}
        void execute(){m_gui.inform("test button pressed");}
        gui::GUI& m_gui;
    };

    // now get the filemenu and add a source button to it.
    gui::SubMenu& filemenu = m_guiMgr->menu().file_menu();
    
    filemenu.addButton("test button", new TestCommand(m_guiMgr->gui()));
}

