// $Header$

#ifdef WIN32  // stupid to prevent compilation on unix platforms 

#include "WinGUI.h"
#include "winGUIostream.h"
#include "resource.h"
#define USECONSOLE  // rather than create a window, that doesn't work now
#include <errno.h>
#include <iostream>
#include <strstream>

static std::streambuf* sbout=0;
static std::streambuf* sberr=0;

//////////////////////////////////////////////////////////////////////
//          private streambuf used by WinGUIostream
class logbuf : public std::streambuf
{
public:
    logbuf (HWND owner);
    virtual ~logbuf ();

    virtual int	overflow (int c);
    virtual int	underflow ();

private:
    HWND _win;

};
//////////////////////////////////////////////////////////////////////
//         WINGUIostream implementations
//////////////////////////////////////////////////////////////////////
//                constructor
WinGUIostream::WinGUIostream(const char* name, HWND owner, HINSTANCE hInstance)
//: ostream(_buf = new logbuf(owner))
: basic_ostream<char>(_buf = new logbuf(owner))
{
 // direct cout and cerr to the console window
    std::streambuf* sbout = std::cout.rdbuf(); 
    std::streambuf* sberr = std::cerr.rdbuf(); 
    std::cout.rdbuf(_buf);
    std::cerr.rdbuf(_buf);

#ifdef USECONSOLE
    AllocConsole();
#else
    _win =  ::CreateDialog( hInstance,MAKEINTRESOURCE(IDD_output), owner, NULL);
#endif
 }
//------------------------------------------------------------------------
WinGUIostream::~WinGUIostream()
{
    std::cout.rdbuf(sbout);
    std::cerr.rdbuf(sberr);
    delete _buf;
}
//------------------------------------------------------------------------
//                           clear

void WinGUIostream::clear()
{
    //TODO: figure out how to clear console
   flush();

}
//////////////////////////////////////////////////////////////////////////////
//            logbuf implementations
//////////////////////////////////////////////////////////////////////////////

const int LOGBUF_SIZE = 32;

//#doc		Initializes the log buffer.
logbuf::logbuf (HWND owner): _win(owner)
{
	setp (0, 0);
	setg (0, 0, 0);
}

logbuf::~logbuf ()
{ }

int logbuf::underflow () { return traits_type::eof (); }

int logbuf::overflow (int c)
{
    if (c != traits_type::eof ())
    {
	char string [3];
	int length = 1;
	string [0] = c;
	string [1] = '\0';
	string [2] = '\0';
	
	if (c == '\n') {
	    string [0] = '\r';
	    string [1] = '\n';
	    ++length;
	}
#ifndef USECONSOLE
	::SendMessage(GetDlgItem(_win,IDC_EDIT1), WM_SETTEXT,0,(LPARAM)string);
	::SetDlgItemText(_win, IDC_EDIT1,string);
	::ShowWindow(_win, SW_RESTORE);
#else
	unsigned long numwritten;
	::WriteConsole(GetStdHandle(STD_OUTPUT_HANDLE), string, length, &numwritten, NULL);
#endif
    }
    return 0;
}


#endif // NT_MSVCPP
