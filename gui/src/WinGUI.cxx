//   $Header$
// Implementation of the GUI interface for MS Windows95/NT
#ifdef WIN32

#include "WinGUI.h"
#include "WinGUIostream.h"
#include "WinScene.h"
#include "resource.h"
#include <iostream>
#include <strstream>
#include "gui/SubMenu.h"

extern "C" { int main(int argn, char *argc[] );} 

static char *szAppName = "WinGUI";
static bool quitting = false;

HINSTANCE   WinGUI::s_instance =0;
WinGUI*     WinGUI::s_gui=0;
HWND        WinGUI::s_hwnd=0;     // windows handle to the window
WinScene*   WinGUI::s_graphics_window=0;
WinScene*   WinGUI::s_2d_window=0;

WinGUIostream* text_window=0;
 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// This enforces the global singleton
namespace gui {
GUI* GUI::createGUI(const char* nameOfApp, const char* title)
{
    if (WinGUI::s_gui!=0 && title!=0 ) ::SetWindowText(WinGUI::s_hwnd, title);
    if (WinGUI::s_gui==0)  WinGUI::s_gui =new WinGUI(nameOfApp, title);
    return WinGUI::s_gui;
}
GUI* GUI::instance()
{
    return  (WinGUI::s_gui)? WinGUI::s_gui :  createGUI("windows");
}
}
bool GUI::running=false; // needed to check that display enabled 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//              constructor
WinGUI::WinGUI(const char* nameOfApp, const char* title)
:   GUI(title)
,   m_menuId(0)
,   m_text_window(0)
{
    // now create the main window, using the registered windows class
    if( s_hwnd==0 )
        s_hwnd = CreateWindow(szAppName, 
        title,
        (WS_OVERLAPPED | WS_CAPTION | WS_MINIMIZEBOX |
        WS_MAXIMIZEBOX | WS_THICKFRAME | WS_SYSMENU),//WS_OVERLAPPEDWINDOW,
        20,  20, //CW_USEDEFAULT, 0,
        500, 500, //CW_USEDEFAULT, 0,
        NULL, NULL, 
        s_instance, 
        NULL);
    else setTitle(title);
    
    // create toplevel menu bar, set handles to it
    m_menuBar = m_currentMenu = ::CreateMenu ();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
WinGUI::~WinGUI()
{
    CommandList::iterator pcmd = m_owned_commands.begin();
    while (pcmd != m_owned_commands.end() )
        delete (*pcmd++);
    delete m_text_window;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
WinGUI::Menu*
WinGUI::beginPullDownMenu(const char * name, WinGUI::Menu* parent)
{
    m_currentMenu = ::CreatePopupMenu();
    ::AppendMenu( parent? (HMENU)parent: m_menuBar, MF_POPUP, (UINT)m_currentMenu, name);
    return (Menu*)m_currentMenu;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void
WinGUI::restorePullDownMenu(Menu * oldMenu)
{
    m_currentMenu = (HMENU)oldMenu;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void
WinGUI::endPullDownMenu(){m_currentMenu = m_menuBar;}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void
WinGUI::addToMenu(const char* name, Command* cmd)
{
    ::AppendMenu(m_currentMenu, MF_STRING, cmdOffset+(m_menuId++), name);
    m_commands.push_back(cmd);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class WinGUI::Toggle : public GUI::Toggle
{ 
    friend class WinGUI;
    Toggle(HMENU hMenu, UINT uPos,  Command* cmd1, Command*cmd2, bool state)
        : m_hMenu(hMenu), m_uPos(uPos)
        , m_cmd1(cmd1)
        , m_cmd2(cmd2)
        , m_state(state)
    {}
    void execute() {
        if( m_state){  m_cmd2->execute();   m_state=FALSE;}
        else	    {  m_cmd1->execute();   m_state=TRUE;	}
        ::CheckMenuItem(m_hMenu, m_uPos, m_state? MF_CHECKED : MF_UNCHECKED);
    }
    bool state()const{return m_state;}

    HMENU m_hMenu; UINT m_uPos;
    Command *m_cmd1;
    Command *m_cmd2;
    bool m_state;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void WinGUI::menuSeparator()
{
    ::AppendMenu(m_currentMenu, MF_SEPARATOR, 0, 0);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GUI::Toggle* WinGUI::addToggleToMenu(const char* name, bool state,  Command* cmd1, Command* cmd2)
{
    UINT uPos = cmdOffset+(m_menuId++);
    GUI::Toggle* t = new Toggle(m_currentMenu, uPos, cmd1,cmd2,state);
    m_commands.push_back(t);
    m_owned_commands.push_back(t);
    ::AppendMenu(m_currentMenu, MF_STRING | (state? MF_CHECKED : MF_UNCHECKED ), uPos,  name);
    return t;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void WinGUI::execute(WPARAM cmdindex)
{
    m_commands[cmdindex-cmdOffset]->execute();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void WinGUI::setTitle(const char* newtitle)
{
    ::SetWindowText(s_hwnd, newtitle);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                     commands
void WinGUI::start()
{
    static bool first=true;
    if( first ) {
        // first-time setup. Allow message loop to be recalled as a result of a message
        ::SetMenu(s_hwnd,m_menuBar);
    
        int cmdShow=SW_SHOWDEFAULT	;
        ::ShowWindow(s_hwnd, cmdShow);
        first = false;
    }
    running = true;

}
void CALLBACK MyTimerProc(HWND hwnd, UINT uMsg, UINT idEvnet, DWORD dwTIme)
{
	GUI::running = false;
}
void WinGUI::run(int pause_interval)
{
    start(); // make sure windows set up
    MSG		msg;
	const UINT uIDEevnt=1;
	if(pause_interval==0) return; 
	if( pause_interval>0 ) 	::SetTimer(s_hwnd, uIDEevnt, pause_interval, MyTimerProc);
	
    // note that the loop can be stopped, with return from this function,
    // by setting GUI::running to false.
    while( running && ::GetMessage(&msg, NULL, 0, 0)) {
        ::TranslateMessage(&msg);
        ::DispatchMessage(&msg);
    }
	if( pause_interval>0) ::KillTimer(s_hwnd, uIDEevnt);
    // unless it is a quit message, reset when leave so main loop does not exit 
    if( !quitting ) running = true;  
    
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void WinGUI::processMessages()
{
    // used to check that messages are processed
    MSG msg;
    
    while( ::PeekMessage(&msg, 0,0,0,PM_REMOVE))
    {
        ::TranslateMessage(&msg);
        ::DispatchMessage(&msg);
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void WinGUI::quit(){::PostQuitMessage(0); quitting=true; running=false;}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                     graphics
// return pointer to a SceneControl sub class that will write into the main window
gui::SceneControl*
WinGUI::graphicsWindow(float size, int initial_view)
{
    if( s_graphics_window==0)  
        return s_graphics_window = new WinScene(s_hwnd, size, initial_view);

    // second window: use same class
    HWND hwnd = ::CreateWindow(szAppName, 
        "plots",
        ( WS_OVERLAPPED | WS_CAPTION | WS_MINIMIZEBOX | WS_POPUP | WS_SYSMENU
        | WS_MAXIMIZEBOX | WS_THICKFRAME | WS_VISIBLE),
        50,  50, //CW_USEDEFAULT, 0,
        500, 500, //CW_USEDEFAULT, 0,
        NULL, NULL, 
        s_instance, 
        NULL);
    return s_2d_window = new WinScene(hwnd, size);

}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//          text output
std::ostream*
WinGUI::textWindow(const char* name)
{
    if( m_text_window==0) {
        if( text_window ==0 ) text_window = new WinGUIostream( name, s_hwnd, s_instance);
        m_text_window = text_window;
    }
    return m_text_window;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//         file name dialog
#include <commdlg.h>
char *
WinGUI::askForFileName(const char* startDir , const char* label, const char* filter)
{
    static OPENFILENAME ofn;

    char szDirName[256];
    char szFile[256], szFileTitle[256];
    char  szFilter[256];
    
    strncpy(szFilter, filter, sizeof(szFilter));
    szFile[0] = '\0';
    
    // get startDir, strip after last "\" character, since it is the last file;
    strncpy(szDirName, startDir, sizeof(szDirName));
    int lpos=0; 
    for( int pos = 0; szDirName[pos]!='\0' ; ++pos ) {
        if( szDirName[pos]=='\\')lpos=pos;
    }
    if( lpos) szDirName[lpos]='\0';
    
    UINT cbString = strlen(szFilter);
    char  chReplace = szFilter[cbString - 1]; /* retrieve wild character */
    
    for (UINT i = 0; szFilter[i] != '\0'; i++) {
        if (szFilter[i] == chReplace)
            szFilter[i] = '\0';
    }
    
    
    /* Set all structure members to zero. */
    
    memset(&ofn, 0, sizeof(OPENFILENAME));
    
    ofn.lStructSize = sizeof(OPENFILENAME);
    ofn.hwndOwner = s_hwnd;
    ofn.lpstrFilter = szFilter;
    ofn.nFilterIndex = 1;
    
    ofn.lpstrTitle= label; 
    ofn.lpstrFile= szFile;
    ofn.nMaxFile = sizeof(szFile);
    ofn.lpstrFileTitle = szFileTitle;
    ofn.nMaxFileTitle = sizeof(szFileTitle);
    ofn.lpstrInitialDir = szDirName;
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_HIDEREADONLY ; //OFN_SHOWHELP | OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;
    
    int result =::GetOpenFileName(&ofn);
    
    // convert backslashes
    /*   for (i = 0; ofn.lpstrFile[i] != '\0'; i++) {
    if (ofn.lpstrFile[i] == '\\')
    ofn.lpstrFile[i] = '/';
}*/
    return result ? ofn.lpstrFile : 0;
    
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                   askuser dialog box
static char userAnswer[80];
static const char* defString;
static const char* prmpt;
BOOL APIENTRY AskUserProc(HWND hdlg, UINT wMsg,
                          WPARAM wParam, LPARAM lParam)
{
    HWND textBox = ::GetDlgItem(hdlg,IDC_ASKUSER);
    switch( wMsg )
    {
    case WM_INITDIALOG:
        // set prompt in title bar and default in text box
        ::SendMessage(textBox, WM_SETTEXT, 0, (LPARAM)defString);
        
        // this should select the text just set
        ::SendMessage(textBox,EM_SETSEL,0, -1);
        
        ::SendMessage(hdlg,WM_SETTEXT, 0, (LPARAM)prmpt);
        return TRUE;
        
    case WM_COMMAND:
        switch (wParam )
        {
        case IDOK:
            ::SendMessage(textBox, WM_GETTEXT, sizeof(userAnswer), (LPARAM)userAnswer);
            ::EndDialog(hdlg, TRUE);
            return TRUE;
        case IDCANCEL:
            ::EndDialog(hdlg, FALSE);
            return TRUE;
        }
        return TRUE;
        default:
            return FALSE;
    }
}
char*
WinGUI::askUser(const char* promptString, const char* defaultString)
{
    userAnswer[0]='\0';
    defString = defaultString;
    prmpt = promptString;
    int status = ::DialogBox( s_instance,MAKEINTRESOURCE(IDD_ASKUSER), s_hwnd, AskUserProc);
    return  userAnswer;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                   message
void
WinGUI::inform(const char* msg)
{
    ::MessageBox(s_hwnd, (char*)msg, "Inform", MB_OK);
}

//========================================================================
//               Windows routines
//========================================================================

LRESULT CALLBACK   WndProc (HWND hWnd, unsigned Message,
                          WPARAM wParam, LONG lParam)
{
    // decide which window this is, set appropriate scene to message
    WinScene* display = ( hWnd==WinGUI::s_hwnd  || WinGUI::s_hwnd==0)?
	WinGUI::s_graphics_window :
	WinGUI::s_2d_window;
    switch(Message)
    {
        
    case WM_COMMAND:
        WinGUI::s_gui->execute(wParam);
        break;
        
    case WM_PAINT:
        {
            PAINTSTRUCT ps;
            HDC hdc=  BeginPaint( hWnd, &ps);
            if(display)
                display->redisplay(hdc);
            EndPaint(hWnd, &ps);
        }
        break;
    case WM_CHAR:
        {	// single non-special key
            char c = static_cast<char>(wParam);
            
            // check if is is registered with the Menu
            if( Menu::instance()->strike(c)) break;
            
            // no, just pass it to the display
            char cmd[2] = {c,0};
            if( display->do_command(cmd)>0 )
                display->redisplay();
        }
        break;
    case WM_KEYDOWN:
        { char cmd[2] = {0,0};
        switch (wParam){
        case VK_LEFT: cmd[1]=75; break;
        case VK_RIGHT:cmd[1]=77; break;
        case VK_HOME: cmd[1]=71; break;
        case VK_UP:   cmd[1]=72; break;
        case VK_DOWN: cmd[1]=80; break;
        case VK_PRIOR:cmd[1]=73; break;
        case VK_NEXT: cmd[1]=81; break;
        }
        if (!cmd[1] )break;
        if( display->do_command(cmd)>0 )
            display->redisplay();
        
        }
        break;
    case WM_LBUTTONDOWN:
        {
            int xPos = LOWORD(lParam);  // horizontal position of cursor
            int	yPos = HIWORD(lParam);  // vertical position of cursor
            display->mouseDown(xPos,yPos);
        }
        break;
    case WM_LBUTTONUP:
        {
            int xPos = LOWORD(lParam);  // horizontal position of cursor
            int	yPos = HIWORD(lParam);  // vertical position of cursor
            display->mouseUp(xPos,yPos);
        }
        break;
    case WM_RBUTTONUP :
        {

	    RECT winpos;
	    ::GetWindowRect(hWnd, &winpos);
	    ::TrackPopupMenu(
		(HMENU)(display->menu()->tag()),  // handle to shortcut menu
		TPM_LEFTALIGN | TPM_TOPALIGN,   // screen-position and mouse-button flags
		winpos.left+LOWORD(lParam),       // horizontal position, in screen coordinates
		winpos.top+HIWORD(lParam),        // vertical position, in screen coordinates
		0,		// reserved, must be zero
		hWnd,           // handle to owner window
		0);
        }
        break;
    case WM_ACTIVATEAPP:
        {
            int fActive = LOWORD(wParam);           // activation flag 
            if( fActive)
                ::SetCursor(::LoadCursor(NULL, IDC_ARROW));
            break;
        }
    case WM_DESTROY:
        ::PostQuitMessage(0);
        quitting=true;
        GUI::running=false;
        break;
    case WM_SETCURSOR:  // mouse comes over window--make it the normal arrow
        {
            ::SetCursor(::LoadCursor(NULL, IDC_ARROW));
            break;
        }
    default:
        return ::DefWindowProc(hWnd, Message, wParam, lParam);
    }
    return 0;
}  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int WINAPI WinMain(HINSTANCE hThisInstance, HINSTANCE hPrevInstance,
                   LPSTR lpszCmdLine, int nCmdShow)
{
    WinGUI::s_instance =  hThisInstance;  //need the instance

    // declare, and create a windows class for the main window
    WNDCLASS	wndClass;

    if (!hPrevInstance)    {
        wndClass.style		= CS_HREDRAW | CS_VREDRAW ;
        wndClass.lpfnWndProc    = WndProc;
        wndClass.cbClsExtra	= 0;
        wndClass.cbWndExtra	= 0;
        wndClass.hInstance	= hThisInstance;
        wndClass.hIcon		= wndClass.hIcon= LoadIcon(hThisInstance,MAKEINTRESOURCE(IDI_ICON1));;
        wndClass.hCursor	= NULL;
        wndClass.hbrBackground	= (HBRUSH)GetStockObject (LTGRAY_BRUSH);
        wndClass.lpszMenuName	= NULL;
        wndClass.lpszClassName  = szAppName;
        
        if( !RegisterClass(&wndClass) ) return FALSE;
    }
    
    // create the main window to be used later
    WinGUI::s_hwnd = CreateWindow(szAppName, 
        "",
        (WS_OVERLAPPED | WS_CAPTION | WS_MINIMIZEBOX |
        WS_MAXIMIZEBOX | WS_THICKFRAME | WS_SYSMENU),//WS_OVERLAPPEDWINDOW,
        20,  20, //CW_USEDEFAULT, 0,
        500, 500, //CW_USEDEFAULT, 0,
        NULL, NULL, 
        WinGUI::s_instance, 
        NULL);

    // setup the console window here
    text_window = new WinGUIostream( "wingui", WinGUI::s_hwnd, WinGUI::s_instance);
    
    //  parse command line by hand
    const unsigned MAX_TOKENS=20;
    char * _argv[MAX_TOKENS];
    _argv[0] = szAppName;
    int argc = 1;
    
    char seps[]   = " \t\n";
    
    /* Establish string and get the first token: */
    char * token = strtok( lpszCmdLine, seps );
    while(token != NULL && argc < MAX_TOKENS)
    {
        /* While there are tokens in "string" */
        _argv[argc++] = token;
        /* Get next token: */
        token = strtok( NULL, seps );
    }
//    ::UserMain(argc, _argv);
    ::main(argc, _argv);
    return 0;
}
#endif
