// Test program for gui graphical display
// $Header$


#if defined(_DEBUG) && defined (_WIN32) //boilerplate for memory leak check
#  include <CRTDBG.H>
#  pragma comment(lib, "msvcrtd.lib")
#endif
#include "Demo.h"
// gui event graphics include files
#include "gui/GuiMgr.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//         demo main 
int main(int /*argc */, char** /*argv*/)
{
  //force into link  
    if(0)gui::GUI::createGUI();


    // get pointer to the GuiMgr instance 
    gui::GuiMgr* guiMgr=gui::GuiMgr::instance();

    // create the demonstration guy.
    Demo demo(guiMgr->display());

    // start the message loop
    guiMgr->menu().run();

    // now the simulated event loop
    do {
        // begin event processing
        guiMgr->begin_event(); 
  
        // do something!
        demo.process();

        // end event processing
        guiMgr->end_event(); 
    }while( !guiMgr->done() );


#if defined(_DEBUG) && defined (_WIN32)
    // this to verify that memory has been released
    // _CrtSetDbgFlag( _CrtSetDbgFlag( _CRTDBG_REPORT_FLAG ) | _CRTDBG_LEAK_CHECK_DF);
#endif
    
    return 0;
}

