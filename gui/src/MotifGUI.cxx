//     $Header$
//  Author: G. Barrand, T. Burnett
// Motif MotifGUI implementation

#ifndef WIN32  // stupid to prevent compilation in windows 

#include "MotifGUI.h"
#include "Xostream.h"
#include "XScene.h"

#include "gui/Command.h"

#include <Xm/Xm.h>
#include <Xm/Form.h>
#include <Xm/RowColumn.h>
#include <Xm/CascadeB.h>
#include <Xm/PushB.h>
#include <Xm/DrawingA.h>
#include <Xm/FileSB.h>
#include <Xm/MessageB.h>
#include <Xm/SelectioB.h>
#include <Xm/Text.h>
#include <Xm/ToggleB.h>
#include <Xm/Separator.h>


#include <strstream>
#include <iostream>
#include <cstdio>
#include <cassert>

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//          File-scope variables (could be class or object variables)
static XtAppContext app;
static Widget top;      //top level widget
static Widget menubar;  // the menu bar
static Widget currentMenu; // either the menu bar or a pull-down menu

static Widget darea=0;
static GUI* theGUI=0;   // pointer to the only allowed instance
static Arg arglist[5];  // use to set args
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    Implement the GUI-creator of the GUI superclass
namespace gui {
GUI* GUI::createGUI(const char* nameOfApp, const char* title)
{
    if( theGUI) return theGUI; // ensures only one instance
    theGUI =new MotifGUI(nameOfApp, title);
    GUI::instance(theGUI);
    GUI::s_instance = theGUI; //(don't know if needed)
    return theGUI:
}
bool GUI::running=false;

}//namespace gui

// forward declaration of  local utilities

static void  XmLabelSetString   (Widget,char*);
Pixel XWidgetStringPixel (Widget,char*);
char* XmStringToString ( XmString a_cps,int a_number);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//              constructor
MotifGUI::MotifGUI(const char* /*nameOfApp*/, const char* title)
:GUI(title), m_scene(0)
{
    static int argc=0;
    static char** argv; // connect to startup?

    //	 create toplevel widget, set app context
    XtSetArg (arglist[0],XtNtitle,(char*)title);
    top =   XtAppInitialize(&app,"MotifoGUI", NULL,0,  &argc,argv,NULL,arglist,1);
    if(!top)   exit(-1);

    // create a form that will hold the menubar and drawing area
    Widget form          = XmCreateForm(top,"form",arglist,0);
    XtManageChild (form);

    // Create the menu bar, which will contain pull-down menus, or simple commands
    menubar  = XmCreateMenuBar(form,"menubar",arglist,0);
    XtManageChild (menubar);

    currentMenu = menubar;

    // Create a drawing area
    XtSetArg      (arglist[0],XmNwidth, 500);
    XtSetArg      (arglist[1],XmNheight,500);
    XtSetArg      (arglist[2],XmNbackground,XWidgetStringPixel (top,"lightgrey"));
    darea         = XmCreateDrawingArea (form,"darea",arglist,3);
    XtManageChild (darea);

    // Glue menubar and drawing area in the form.
    XtSetArg      (arglist[0],XmNtopAttachment,XmATTACH_FORM);
    XtSetArg      (arglist[1],XmNleftAttachment,XmATTACH_FORM);
    XtSetArg      (arglist[2],XmNrightAttachment,XmATTACH_FORM);
    XtSetValues   (menubar,arglist,3);

    XtSetArg      (arglist[0],XmNtopAttachment,XmATTACH_WIDGET);
    XtSetArg      (arglist[1],XmNtopWidget,menubar);
    XtSetArg      (arglist[2],XmNleftAttachment,XmATTACH_FORM);
    XtSetArg      (arglist[3],XmNrightAttachment,XmATTACH_FORM);
    XtSetArg      (arglist[4],XmNbottomAttachment,XmATTACH_FORM);
    XtSetValues   (darea,arglist,5);

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void MotifGUI::setTitle(const char* newtitle)
{
    XtSetArg (arglist[0],XtNtitle,(char*)newtitle);
    XtSetValues (top, arglist,1);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//              Menu creation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static void ButtonCallback(Widget , XtPointer clientData, XtPointer /*callData*/)
{
    Command* cmd = (Command*)clientData;
    cmd->execute();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GUI::Menu*
MotifGUI::beginPullDownMenu(const char * menuName, GUI::Menu* parent)
{
    //if parent is there, use it.
    Widget pmen = parent? (Widget)parent : menubar;
    Widget casc   = XmCreateCascadeButton(pmen,"casc",arglist,0);
    XtManageChild (casc);

    Widget pull   = XmCreatePulldownMenu (pmen,"pull",arglist,0);
    XtSetArg      (arglist[0], "subMenuId", pull);
    XtSetValues   (casc,arglist,1);
    XmLabelSetString (casc,(char*)menuName);

    // allow to be torn off?
    XtSetArg( arglist[0], XmNtearOffModel, XmTEAR_OFF_ENABLED);
    XtSetValues( pull, arglist, 1);

    currentMenu = pull;
    return (GUI::Menu*) pull;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void MotifGUI::restorePullDownMenu(GUI::Menu* m)
{
    currentMenu = (Widget)m;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void MotifGUI::addToMenu(const char* buttonName, Command* cmd)
{

  if( currentMenu == menubar ) {
    // kluge to allow single button menu??? (hard to debug otherwise)
    beginPullDownMenu(buttonName);
    addToMenu(buttonName, cmd);
    endPullDownMenu();
    return;
  }

    Widget button = XmCreatePushButton (currentMenu, "button", arglist,0);

    XtManageChild (button);
    XtAddCallback (button, XmNactivateCallback, ButtonCallback,(XtPointer)cmd);
    XmLabelSetString (button,(char*)buttonName);

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class MotifGUI::Toggle : public GUI::Toggle
// private class handles toggle button
{
public:
    Toggle( Command* cmd1, Command* cmd2 ,bool state, Widget button)
	:  m_set_cmd(cmd1)
	,  m_unset_cmd(cmd2)
	, m_state(state)
        , m_button(button)
    {}
    void execute() {
	if( m_state){
            m_unset_cmd->execute();   m_state=false;
	}
	else	    {
            m_set_cmd->execute();   m_state=true;
	}
	// should just set correct widget state if not called from X 
	XtSetArg(arglist[0], XmNset, m_state );	
	XtSetValues(m_button, arglist, 1);
    }
  bool state()const{return m_state;}
private:
    Command *m_set_cmd, *m_unset_cmd;
    bool m_state;
    Widget m_button;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GUI::Toggle*  MotifGUI::addToggleToMenu(const char* title, bool state, Command* set, Command* unset)
{
    XtSetArg(arglist[0], XmNset, state );
    Widget button =  XmCreateToggleButton(currentMenu, "button", arglist,1);
    XtManageChild (button);
    GUI::Toggle* t = new Toggle(set,unset,state,button);
    XtAddCallback (button, XmNvalueChangedCallback, ButtonCallback,  (XtPointer)t);
    XmLabelSetString (button, (char*)title);
    return t;

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void MotifGUI::menuSeparator()
{
    XtManageChild(
	XmCreateSeparator(currentMenu, "separator", arglist, 0));
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void MotifGUI::endPullDownMenu()
{
    currentMenu = menubar;  // reset, so next command will go directly to menu bar
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                     commands
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static void TimerCallback (  XtPointer a_data, XtIntervalId* timer )
{
	GUI::running = false; //
}

void MotifGUI::run(int pause_interval)
{
	static bool first = true;
	if (first) {
		XtRealizeWidget(top);
		first = false;
	}
	if( pause_interval ==0 ) return; // no loop
	XtIntervalId timer=0;
	XtPointer client_data=0;
  
	if( pause_interval>0) { // set up a timer
		timer = XtAppAddTimeOut(app, pause_interval, &TimerCallback, client_data);
	}
	running = true;
	// the following is equivalent to    XtAppMainLoop(app);
	// execpt that it can exit
    while( running){
		//processMessages();
#if 0		
		XEvent event;
		XtAppNextEvent(app, &event);
		XtDispatchEvent(&event);
#endif
		XtAppProcessEvent (app, XtIMAll);
    }
	if( timer!=0) XtRemoveTimeOut(timer); // in case resume by hand
    running = true; //if this is a secondary loop, so that the main loop doesn't exit 
    
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void MotifGUI::quit()
{
    running = false;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void MotifGUI::processMessages()
{
    while( XtAppPending(app) & XtIMXEvent) {
	XEvent event;
	XtAppNextEvent(app, &event);
	XtDispatchEvent(&event);
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                     graphics
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SceneControl* MotifGUI::graphicsWindow(float size, int initial_view)
{
    Widget drawWidget = darea;  // see what is happening
    m_scene = new XScene(size, initial_view, drawWidget);

    return  m_scene;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                    text output
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
std::ostream* MotifGUI::textWindow(const char* name)
{
    return  new Xostream(  name ) ;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                    file dialog
static int messageAck;
static Widget fileDialog=0;
char* chosenFile; // communicate file chooser result
static void OkCallback (  Widget w ,XtPointer /*a_tag*/, XtPointer a_data )
{
    XmFileSelectionBoxCallbackStruct *data
                = (XmFileSelectionBoxCallbackStruct*)a_data;
    chosenFile = XmStringToString (data->value,0);
    XtUnmanageChild(w);
    messageAck=1;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
char * MotifGUI::askForFileName(const char* /*startDir*/ , const char* label,
				const char* /*boxLabel*/)
{
    if( !fileDialog ){
	// Create a File chooser
	XtSetArg      (arglist[0],XmNautoUnmanage,True);
	fileDialog   = XmCreateFileSelectionDialog (top,"fileDialog",arglist,1);
	XtAddCallback (fileDialog,XmNokCallback  ,OkCallback,(XtPointer)0);
    }
    XtSetArg (arglist[0],XtNtitle,(char*)label);
    XtSetValues (fileDialog, arglist, 1);

    messageAck=0;
    XtManageChild(fileDialog);
//   XtAddGrab(fileDialog, True, False);
    while( !messageAck)   {
	XEvent event;
	XtAppNextEvent(app, &event);
	XtDispatchEvent(&event);
    }
    //?? XtRemoveGrab(fileDialog);

    return chosenFile;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                            message dialog
static Widget messageDialog=0;
static void MessageCB(Widget w, XtPointer /* client_data */,
		      XmAnyCallbackStruct * /*call_data*/)
{
    XtUnmanageChild(w);
    messageAck=1;
}
void
MotifGUI::inform(const char* msg)
{
    if( !messageDialog )   {
	messageDialog=XmCreateMessageDialog(top,"dialog",arglist,0);
	XtAddCallback(messageDialog,XmNokCallback,(XtCallbackProc)MessageCB,NULL);
	XtUnmanageChild(XmMessageBoxGetChild(messageDialog,XmDIALOG_CANCEL_BUTTON));
	XtUnmanageChild(XmMessageBoxGetChild(messageDialog,XmDIALOG_HELP_BUTTON));
    }
    XtSetArg(arglist[0], XmNmessageString,
        XmStringCreateLtoR((char*)msg,XmSTRING_DEFAULT_CHARSET));
    XtSetValues   (messageDialog,arglist,1);
    messageAck=0;
    XtManageChild(messageDialog);
 //  XtAddGrab(messageDialog, True, False);
    while( !messageAck )   {
	XEvent event;
	XtAppNextEvent(app, &event);
	XtDispatchEvent(&event);
    }
    //?? XtRemoveGrab(messageDialog);

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                    prompt dialog
static Widget promptDialog=0;
static char* answer;

static void PromptCB(Widget w, XtPointer /*client_data */,
           XmSelectionBoxCallbackStruct  *call_data)
{
    if( call_data->reason==XmCR_OK )
	XmStringGetLtoR(call_data->value,XmSTRING_DEFAULT_CHARSET,&answer);
    messageAck=1;
    XtUnmanageChild(w);
}

char*
MotifGUI::askUser(const char* promptString,const char* defaultString)
{
    if( !promptDialog )   {
	promptDialog=XmCreatePromptDialog(top,"dialog",arglist,0);

	XtAddCallback(promptDialog,XmNokCallback,(XtCallbackProc)PromptCB,NULL);
	XtAddCallback(promptDialog,XmNcancelCallback,(XtCallbackProc)PromptCB,NULL);

	//doesn't work? XtUnmanageChild(XmMessageBoxGetChild(promptDialog,XmDIALOG_HELP_BUTTON));

    }
    XtSetArg(arglist[0], XmNselectionLabelString,
        XmStringCreateLtoR((char*)promptString,XmSTRING_DEFAULT_CHARSET));
    XtSetArg(arglist[1], XmNtextString,
        XmStringCreateLtoR((char*)defaultString,XmSTRING_DEFAULT_CHARSET));
    XtSetValues   (promptDialog,arglist,2);
    //  XtAddGrab(promptDialog, True, False);
    XtManageChild(promptDialog);
    messageAck=0;
    while( !messageAck )   {
	XEvent event;
	XtAppNextEvent(app, &event);
	XtDispatchEvent(&event);
    }
    //?? XtRemoveGrab(promptDialog);

    return answer;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static void XmLabelSetString (
 Widget This
, char* a_string
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  Arg       arglist[1];
  XmString  cps;
/*.........................................................................*/
  if(!This)     return;
  if(!a_string) return;
  cps          = XmStringLtoRCreate(a_string,XmSTRING_DEFAULT_CHARSET);
  XtSetArg     (arglist[0],XmNlabelString,cps);
  XtSetValues  (This,arglist,1);
  XmStringFree (cps);
}
/***************************************************************************/
Pixel XWidgetStringPixel (
 Widget This
,char* a_string
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  XrmValue  from,to;
/*.........................................................................*/
  if(!This)                     return (Pixel)NULL;
  if(!a_string || !(*a_string)) return (Pixel)NULL;
  from.size    = strlen(a_string)+1;
  from.addr    = (caddr_t)a_string;
  XtConvert    (This,XtRString,&from,XtRPixel,&to);
  if(!to.addr) return (Pixel)NULL;
  return       (*((Pixel*)to.addr));
}

char* XmStringToString (
 XmString a_cps
,int a_number
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  char*             string;
  register int      icount;
  XmStringContext   context;
  Boolean           done;
/*.........................................................................*/
  if(!a_cps)           return 0;
  string               = 0;
  XmStringInitContext  (&context,a_cps);
  icount               = 0;
  done                 = False;
  while(!done)
    {  char*             text    = 0;
       XmStringCharSet   charset = 0;
       XmStringDirection direct;
       Boolean           sep;
       if(XmStringGetNextSegment (context,&text,&charset,&direct,&sep))
         {
           XtFree               (charset);
           if(sep)              done = True;
           if(icount==a_number)
             {
               string = text;
               break;
             }
           icount               ++;
           XtFree               (text);
         }
       else
         done = True;
    }
  XmStringFreeContext  (context);
  return               string;
}


#endif // of NT_MSVCPP
