//     $Id$
//  Author: T. Burnett, G. Barrand
// 
// It declares the class XScene, a specific  display
// for X11

#include "gui/SceneControl.h"
#include <X11/Intrinsic.h>
#include <iostream>
class Xdraw;

using namespace gui;

class XScene : public  SceneControl
   {
 public:

   XScene(float size=100, int initial_view=1, Widget=0);
   // constructor. If no widget, it will create an X-display window
   ~XScene();
   // destructor

   int viewit();
   // Start an event or detector display

   virtual void redisplay();
   // overrides SceneControl to make sure Xdraw is ready

   int processEvent(XEvent& myevent);
   // viewit and Motif use
 private:
   void viewitInstructions(std::ostream& =std::cout);
   int getEvent();
   // used for viewit interface

   Xdraw*       _drawer;   // associated VwDraw object
   Display*     _display; // X display, window
   Window       _window;

   Widget       _widget; // associated Widget if set by client

 };







