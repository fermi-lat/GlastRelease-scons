#ifdef WIN32 // only used for the windows app
// $Header$

#include "WinGui.h"
#include "WinScene.h"
#include "WinGUIostream.h"
#include "gui/Menu.h"
#include "gui/SubMenu.h"
#include "resource.h"

extern "C" { int main(int argn, char *argc[] );} 

static char *szAppName = "WinGUI";


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
        WinGUI::s_quitting=true;
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
    WinGUI::s_text_window = new WinGUIostream( "wingui", WinGUI::s_hwnd, WinGUI::s_instance);
    
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