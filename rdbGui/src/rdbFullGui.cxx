// $Header$
/**
   @file RDBgui.cxx
   Main for RDBgui application

*/
#include "RdbGuiWindow.h"

// Here we begin
int main(int argc,char *argv[])
{

  // Make application
  FXApp application("RdbGUI","Calibration");

  // Open the display
  application.init(argc,argv);

  // Make window
  new RdbGUIWindow(&application);

  // Create the application's windows
  application.create();

  // Run the application
  return application.run();
}

