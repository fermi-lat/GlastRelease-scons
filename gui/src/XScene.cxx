//  $Header$
//   Author: T. Burnett
//
// Implement XScene functions

#ifndef WIN32  // stupid to prevent compilation in windows 
#include <string>

#include "XScene.h"
#include "Xdraw.h"
#include "MotifGUI.h"
#include "gui/Menu.h"  // for key binding
#include "gui/Projector.h"

#include <X11/keysym.h>
#include <assert.h>
#include <iostream>

typedef struct {int reason;  XEvent *event;
        } XmAnyCallbackStruct;

static void ExposeCallback       (Widget,XtPointer,XmAnyCallbackStruct*);
static void InputCallback        (Widget,XtPointer,XmAnyCallbackStruct*);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                   constructor
XScene::XScene (float size, int initial_view, Widget drawWidget)
: SceneControl ( new Projector(_drawer=new Xdraw), size, initial_view )
, _widget(drawWidget)
{

  if(_widget)
  {
      XtAddCallback (_widget,"exposeCallback" ,(XtCallbackProc)ExposeCallback,(XtPointer)this);
      XtAddCallback (_widget,"inputCallback" ,(XtCallbackProc)InputCallback,(XtPointer)this);
  }
  else  {
      // called without a widget--ceate a display and window

      _display = XOpenDisplay ("");       // connect to the server

      int screen = DefaultScreen (_display);  // get a screen
      _window = XCreateSimpleWindow (_display,
                                DefaultRootWindow (_display),
                                50,50,500,500, 2,
                                BlackPixel(_display, screen),
                                WhitePixel(_display, screen)
                                );


      // tell Window Manager what we are doing
      static int argc=0;
      static char **argv;
      XSizeHints myhint = {  PPosition | PSize, 50,50, 500,500};

      XSetStandardProperties (_display, _window,
                           "gui View", "gui",
                           None, argv, argc, &myhint);

      XSelectInput(_display, _window, ButtonPressMask | ButtonReleaseMask | PointerMotionMask
                | KeyPressMask | ExposureMask);

      XMapRaised (_display, _window);      // pop this window up on the screen

      /* this doesnt' work? (THB)
      // Rebind the Keysyms for the arrow keys and some others, Jones 371
      KeySym modlist[1]; // array of modifier Keysyms
#define REBIND(from,to)\
        XRebindKeysym( _display, from, modlist, 0, (const unsigned char *)to,1)
   REBIND( XK_Left,  "b");
   REBIND( XK_Right, "f");
   REBIND( XK_Up,    "p");
   REBIND( XK_Down,  "n");

   REBIND( XK_Home,  "r");
   REBIND( XK_Prior, "s");
   REBIND( XK_Next,  "m");
   REBIND( XK_Begin, "o");
   REBIND( XK_End,   "i");
#undef REBIND
(THB)*/ 
   _drawer->setDisplay(_display,_window);

  }
}
XScene::~XScene()
{
   delete _drawer;
   if(!_widget )
   {  // created these things, must close them down
     if (_window)
      XDestroyWindow (_display, _window);
     if (_display)
      XCloseDisplay (_display);
   }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                   callbacks


static void ExposeCallback ( Widget , XtPointer disp ,XmAnyCallbackStruct*)
{
  ((XScene*)disp)->redisplay();
}
static void InputCallback ( Widget , XtPointer xptr ,XmAnyCallbackStruct* cbInfo)
{
   XScene& disp = *(XScene*)xptr;
   if( !cbInfo->event ) return;
   int ret = disp.processEvent(*(cbInfo->event));
   if( ret !=0 )  {
       // check if is is registered with the Menu
       if( Menu::instance()->strike((char) ret)) return;

       char command[] = {(char)ret,'\0'};
       ret =disp.do_command(command);
       if( ret ==1)
            disp.displaycurrent(1);
       if( ret==2 )
           disp.redisplay();
   }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//               display -- override superclass
void XScene::redisplay()
{
  if( !MotifGUI::running ) return;
  if( _widget )
      _drawer->setDisplay(XtDisplay(_widget),XtWindow(_widget));
  _drawer->resize();
  _drawer->set_defaults();
  SceneControl::redisplay();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//         This is a displayer for X_Windows Graphics
int XScene::getEvent()
{
   static XEvent myevent;
   XNextEvent(_display, &myevent);
   return processEvent(myevent);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// process events in the window that change or require updates
int XScene::processEvent(XEvent& myevent)
{
   char text[2]={'\0','\0'};
   int i = 0;
   static int buttonState = 0;
   static int selectionCleared = 1;
   static Xdraw::Rectangle selection;

   switch (myevent.type)
      {
    case ButtonPress:
      if( myevent.xbutton.button!=1 ) break;
      selection.x = _drawer->xnorm(myevent.xbutton.x);
      selection.y = _drawer->ynorm(myevent.xbutton.y);
      buttonState = 1;
      break;

    case ButtonRelease:      
      if( myevent.xbutton.button!=1 ) break;
      selection.dx = _drawer->xnorm(myevent.xbutton.x) - selection.x;
      selection.dy = _drawer->ynorm(myevent.xbutton.y) - selection.y;
      setViewPort (selection.x, selection.y,
                              selection.x + selection.dx,
                              selection.y + selection.dy);

     text[0] = 'u';    //force update
      i = 1;
      buttonState = 0;
      selectionCleared = 1;
      break;

    case MotionNotify:
      if ( buttonState )
         {

         // set XOR-mode
         _drawer->setXorMode( 1 );

         // delete old rectangle by drawing it again
         if (!selectionCleared)
            _drawer->drawRect( selection );

         // calc new selection
         selection.dx = _drawer->xnorm(myevent.xbutton.x) - selection.x;
         selection.dy = _drawer->ynorm(myevent.xbutton.y) - selection.y;

         // draw new rectangle
         _drawer->drawRect( selection );
         selectionCleared = 0;

         // clear XOR-mode
         _drawer->setXorMode( 0 );
         }
      break;

    case KeyPress:
      KeySym mykey;
      i = XLookupString ((XKeyEvent*)&myevent, text, 10, &mykey, 0);
      if(i==0){ // special key?
	i = 1;
	switch(mykey){ // if one of these, make translation
	case XK_Left:  text[0]='b'; break;
	case XK_Right: text[0]='f'; break;
	case XK_Up:    text[0]='p'; break;
	case XK_Down:  text[0]='n'; break;
	case XK_Home:  text[0]='r'; break;
	case XK_Prior: text[0]='m'; break;
	case XK_Next:  text[0]='s'; break;
	default: i=0;
	}
      }
      break;

    case VisibilityNotify:
      break;

    case Expose:
      {
      text[0]='u'; i=1;
      }
      break;

    case MappingNotify:
       XRefreshKeyboardMapping ((XMappingEvent*)&myevent);
      break;
      }   // end of    switch (myevent.type)


   if ( i == 0 )
      return 0;
   else
      return text[0];

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                   interactive viewer loop
int XScene::viewit()
{
   print_commands(std::cout);
   std::cout << "q to quit\n";
   redisplay();

   for (;;)
   {
     int command =  getEvent();

     switch(command)
     {
       case 0:   break;

       case 'q': return 0;
       case 't': return 1;
       case 'g': return 2;
       default:    //not recognized here: let SceneControl handle it
         char cmdstring[]={command,'\0'};
         int j= do_command(cmdstring);
         if(j==1) displaycurrent(1);
         if(j==2) redisplay();
      }        // end of switch (command)

   }   // end of the for (;;)

   return 0;
}

#endif // NT_MSVCPP
